# ======================================================================
# == GEO/SRA download rules ==
# ======================================================================
# Optional pre-processing module: given a GSE accession, downloads all
# associated FASTQs from SRA and produces a samples_csv-compatible
# samplesheet (sample_name, sample_run, sample_path).
#
# This module is NOT wired into the main pipeline DAG: df_sample is parsed
# from config["samples_csv"] at workflow parse-time, before any rule runs,
# so it cannot depend on an output that this module produces later.
# Usage is therefore a two-step process:
#   1. snakemake geo_samplesheet   (fetches FASTQs, writes samplesheet.csv)
#   2. point config["samples_csv"] at the generated samplesheet.csv, then
#      run the pipeline as usual (e.g. cellranger_multirun_default).
#
# All outputs are written under the main pipeline's config["out_location"]
# (i.e. inside the results folder), namespaced under "geo_download/" to avoid
# colliding with the rest of the pipeline's outputs. Envs are namespaced
# under config["geo_download"].

GEO = config["geo_download"]
GEO_OUT = config["out_location"] + "geo_download/"
SRA_CACHE = GEO_OUT + "sra_cache/"

# ---------------------------------------------------------------------------- #
checkpoint fetch_metadata:
    '''
    Query NCBI SRA via pysradb to retrieve all SRR run accessions for a GSE.
    Outputs:
      - srr_list.txt  : one SRR accession per line (used to drive downstream rules)
      - metadata.tsv  : full run-level metadata table
    '''
    output:
        srr_list = GEO_OUT + "{gse}/metadata/srr_list.txt",
        metadata  = GEO_OUT + "{gse}/metadata/metadata.tsv"
    params:
        gse = lambda wildcards: wildcards.gse
    conda:
        GEO["env_metadata"]
    log:
        "logs/fetch_metadata/{gse}.log"
    benchmark:
        "benchmarks/fetch_metadata/{gse}.tsv"
    script:
        "../scripts/fetch_srr_from_gse.py"


# ---------------------------------------------------------------------------- #
rule prefetch:
    '''
    Download SRA file for a given SRR accession.
    Uses NCBI prefetch which supports resume and integrity checking.
    Output: {sra_cache}/{srr}/{srr}.sra
    '''
    output:
        sra = SRA_CACHE + "{srr}/{srr}.sra"
    params:
        cache_dir = SRA_CACHE
    conda:
        GEO["env_sratools"]
    log:
        "logs/prefetch/{srr}.log"
    benchmark:
        "benchmarks/prefetch/{srr}.tsv"
    shell:
        '''
        prefetch {wildcards.srr} \
            --output-directory {params.cache_dir} \
            --max-size 200g \
            2> {log}
        '''


# ---------------------------------------------------------------------------- #
rule fasterq_dump:
    '''
    Convert a downloaded SRA file to FASTQ using fasterq-dump.
    --include-technical is required for 10x data: the barcode+UMI (R1) is stored as a
    technical read in SRA and silently skipped without this flag, leaving only the cDNA (R2).
    --split-files writes one file per read position: _1, _2, ... _N.
    Outputs a sentinel .done file; actual FASTQs are in the same directory.
    '''
    input:
        sra = rules.prefetch.output.sra
    output:
        done = touch(GEO_OUT + "{gse}/fastq/{srr}/.fasterq.done")
    params:
        outdir = GEO_OUT + "{gse}/fastq/{srr}",
        core = config["set-threads"]["fasterq_dump"]
    conda:
        GEO["env_sratools"]
    log:
        "logs/fasterq_dump/{gse}/{srr}.log"
    benchmark:
        "benchmarks/fasterq_dump/{gse}/{srr}.tsv"
    shell:
        '''
        mkdir -p {params.outdir}
        fasterq-dump {input.sra} \
            --split-files \
            --include-technical \
            --outdir {params.outdir} \
            --threads {params.core} \
            2> {log}
        '''

# ---------------------------------------------------------------------------- #
rule compress_fastq:
    '''
    Compress all FASTQ files in the SRR directory using pigz (parallel gzip).
    Runs after fasterq_dump sentinel is written.
    '''
    input:
        done = rules.fasterq_dump.output.done
    output:
        compressed = touch(GEO_OUT + "{gse}/fastq/{srr}/.compressed.done")
    params:
        outdir = GEO_OUT + "{gse}/fastq/{srr}",
        core = config["set-threads"]["compress_fastq"]
    conda:
        GEO["env_sratools"]
    log:
        "logs/compress_fastq/{gse}/{srr}.log"
    benchmark:
        "benchmarks/compress_fastq/{gse}/{srr}.tsv"
    shell:
        '''
        find {params.outdir} -maxdepth 1 -name "*.fastq" | xargs -r pigz -f -p {params.core} 2> {log}
        '''

# ---------------------------------------------------------------------------- #
rule get_sra_read_layout:
    '''
    Read labels (R1/R2/I1/I2) as declared in the SRA metadata: READ_SPEC, original
    submitted filenames, or sra-stat Technical/Biological. Metadata only, possibly empty;
    verify_read_layout checks it against the reads and writes the final layout.
    '''
    input:
        sra = rules.prefetch.output.sra
    output:
        layout = GEO_OUT + "{gse}/read_layout/{srr}.metadata.json"
    conda:
        GEO["env_sratools"]
    log:
        "logs/get_sra_read_layout/{gse}/{srr}.log"
    benchmark:
        "benchmarks/get_sra_read_layout/{gse}/{srr}.tsv"
    script:
        "../scripts/get_sra_read_layout.py"


# ---------------------------------------------------------------------------- #
rule locate_barcode_whitelists:
    '''
    Find the 10x barcode whitelists of the cellranger the pipeline runs
    (<install>/lib/python/cellranger/barcodes) and write their directory to a file.
    Uses the same binary name as the count/mkref rules (cellranger10): change it
    here too when the pipeline moves to another cellranger version.
    Not a localrule on purpose: a named env such as env_cellranger="demux" is resolved
    with the profile's conda-base-path only in cluster jobs; a local job uses the conda
    snakemake runs from, where the env may not exist.
    '''
    output:
        whitelist_dir = GEO_OUT + "barcode_whitelist_dir.txt"
    conda:
        config["env_cellranger"]
    log:
        "logs/locate_barcode_whitelists.log"
    shell:
        '''
        exec 2> {log}
        command -v cellranger10 >&2 \
            || {{ echo "cellranger10 not on PATH: env_cellranger not activated? (CONDA_PREFIX=${{CONDA_PREFIX:-unset}})" >&2; exit 1; }}
        cr=$(readlink -f "$(command -v cellranger10)")
        dir=$(dirname "$(dirname "$cr")")/lib/python/cellranger/barcodes
        [[ -d "$dir" ]] || {{ echo "No barcode whitelists at $dir (cellranger10 -> $cr)" >&2; exit 1; }}
        echo "cellranger10 -> $cr" >&2
        echo "$dir" > {output.whitelist_dir}
        '''


# ---------------------------------------------------------------------------- #
rule verify_read_layout:
    '''
    Final R1/R2/I1/I2 labels for one SRR. The first reads of each FASTQ are matched
    against the 10x barcode whitelists to find the barcode read (R1), which is used to
    verify the SRA metadata labels or to label the reads when metadata is missing.
    Fails with an explicit message when the layout cannot be established; see
    config["geo_download"]["read_layout_override"] and the script docstring.
    '''
    input:
        metadata_layout = rules.get_sra_read_layout.output.layout,
        compressed      = rules.compress_fastq.output.compressed,
        metadata        = rules.fetch_metadata.output.metadata,
        whitelist_dir   = rules.locate_barcode_whitelists.output.whitelist_dir
    output:
        layout = GEO_OUT + "{gse}/read_layout/{srr}.json"
    params:
        fastq_dir             = GEO_OUT + "{gse}/fastq",
        n_reads               = GEO.get("read_layout_check", {}).get("n_reads", 10000),
        min_barcode_match     = GEO.get("read_layout_check", {}).get("min_barcode_match", 0.5),
        allow_length_fallback = GEO.get("read_layout_check", {}).get("allow_length_fallback", False),
        overrides             = GEO.get("read_layout_override") or {}
    resources:
        mem_mb = 4000
    conda:
        config["env_bioinfo"]
    log:
        "logs/verify_read_layout/{gse}/{srr}.log"
    benchmark:
        "benchmarks/verify_read_layout/{gse}/{srr}.tsv"
    script:
        "../scripts/verify_read_layout.py"


# ---------------------------------------------------------------------------- #
def get_srrs_for_gse(wildcards):
    srr_file = checkpoints.fetch_metadata.get(gse=wildcards.gse).output.srr_list
    with open(srr_file) as f:
        srrs = [line.strip() for line in f if line.strip()]
    return expand(
        rules.compress_fastq.output.compressed,
        gse=wildcards.gse,
        srr=srrs
    )


def get_layouts_for_gse(wildcards):
    srr_file = checkpoints.fetch_metadata.get(gse=wildcards.gse).output.srr_list
    with open(srr_file) as f:
        srrs = [line.strip() for line in f if line.strip()]
    return expand(
        rules.verify_read_layout.output.layout,
        gse=wildcards.gse,
        srr=srrs
    )


rule group_fastqs_by_sample:
    '''
    Symlink compressed FASTQs into per-sample directories with SpaceRanger/CellRanger-compatible naming.
    Groups SRR runs by experiment (SRX = one biological sample); each SRR becomes a lane.
    Read types come from the verified per-SRR layouts (verify_read_layout). Before linking,
    every lane must have all FASTQs labelled and exactly one R1 and one R2; otherwise the
    rule fails, lists every failing SRR and groups nothing.
    Output: {out_location}/{gse}/grouped/{sample_name}/{sample}_S1_L00N_{R1|R2|I1|I2}_001.fastq.gz
    '''
    input:
        metadata   = rules.fetch_metadata.output.metadata,
        compressed = get_srrs_for_gse,
        layouts    = get_layouts_for_gse
    output:
        done = touch(GEO_OUT + "{gse}/grouped/.done")
    params:
        fastq_dir = GEO_OUT + "{gse}/fastq",
        out_dir   = GEO_OUT + "{gse}/grouped"
    conda:
        GEO["env_metadata"]
    log:
        "logs/group_fastqs/{gse}.log"
    benchmark:
        "benchmarks/group_fastqs/{gse}.tsv"
    script:
        "../scripts/group_fastqs.py"


# ---------------------------------------------------------------------------- #
rule make_samplesheet:
    '''
    Generate a samplesheet CSV from the grouped FASTQ directories.
    Columns: sample_name, sample_run, sample_path — directly compatible with
    this pipeline's config["samples_csv"] input format.
    '''
    input:
        done = rules.group_fastqs_by_sample.output.done
    output:
        csv = GEO_OUT + "{gse}/samplesheet.csv"
    params:
        grouped_dir = GEO_OUT + "{gse}/grouped"
    conda:
        GEO["env_metadata"]
    log:
        "logs/make_samplesheet/{gse}.log"
    benchmark:
        "benchmarks/make_samplesheet/{gse}.tsv"
    script:
        "../scripts/make_samplesheet.py"
