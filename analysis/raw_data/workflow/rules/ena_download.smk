# ======================================================================
# == ENA download rules ==
# ======================================================================
# Optional pre-processing module: given an ENA/SRA project accession (e.g. PRJEB39323),
# downloads the FASTQs submitted to ENA and produces a samples_csv-compatible
# samplesheet (sample_name, sample_run, sample_path).
#
# Use it instead of geo_download when the reads are not in SRA (e.g. GSE118257: GEO/SRA
# lists SRP178139 with no runs, the reads live in ENA PRJEB39323). It downloads the
# submitted files as they are (no prefetch/fasterq-dump), so it only supports projects whose
# submitted files are FASTQs named like Cell Ranger's (<prefix>_S1_L001_R1_001.fastq.gz).
#
# Like geo_download it is NOT wired into the main pipeline DAG (samples_csv is parsed at
# workflow parse-time). Two-step usage:
#   1. snakemake ena_samplesheet   (fetches FASTQs, writes samplesheet.csv)
#   2. point config["samples_csv"] at the generated samplesheet.csv, then run the pipeline.
#
# Outputs are written under config["out_location"] + "ena_download/<project>/".

ENA = config["ena_download"]
ENA_OUT = config["out_location"] + "ena_download/"

# ---------------------------------------------------------------------------- #
checkpoint ena_fetch_metadata:
    '''
    Query the ENA portal API for all runs of a project.
    Outputs:
      - run_list.txt : one run accession per line (used to drive downstream rules)
      - metadata.tsv : run-level table (run, experiment, sample, sample_title, submitted_ftp, submitted_md5, ...)
    '''
    output:
        run_list = ENA_OUT + "{project}/metadata/run_list.txt",
        metadata = ENA_OUT + "{project}/metadata/metadata.tsv"
    params:
        project = lambda wildcards: wildcards.project
    wildcard_constraints:
        project = r"[A-Za-z0-9]+"
    conda:
        ENA["env_metadata"]
    log:
        "logs/ena_fetch_metadata/{project}.log"
    benchmark:
        "benchmarks/ena_fetch_metadata/{project}.tsv"
    script:
        "../scripts/fetch_ena_runs.py"


# ---------------------------------------------------------------------------- #
rule ena_download_fastq:
    '''
    Download the submitted FASTQ files of one run over HTTPS (resumable) and check their MD5.
    Output: a sentinel .done file; the FASTQs are in the same directory.
    '''
    input:
        metadata = rules.ena_fetch_metadata.output.metadata
    output:
        done = touch(ENA_OUT + "{project}/fastq/{run}/.done")
    params:
        outdir = ENA_OUT + "{project}/fastq/{run}",
        run    = lambda wildcards: wildcards.run
    wildcard_constraints:
        project = r"[A-Za-z0-9]+",
        run     = r"[A-Za-z]+[0-9]+"
    resources:
        ena_connections = 1   # capped globally by `resources: ena_connections` in the profile
    conda:
        ENA["env_metadata"]
    log:
        "logs/ena_download_fastq/{project}/{run}.log"
    benchmark:
        "benchmarks/ena_download_fastq/{project}/{run}.tsv"
    script:
        "../scripts/download_ena_fastq.py"


# ---------------------------------------------------------------------------- #
def get_ena_downloads(wildcards):
    run_file = checkpoints.ena_fetch_metadata.get(project=wildcards.project).output.run_list
    with open(run_file) as f:
        runs = [line.strip() for line in f if line.strip()]
    return expand(
        rules.ena_download_fastq.output.done,
        project=wildcards.project,
        run=runs
    )


rule ena_group_fastqs:
    '''
    Symlink the downloaded FASTQs into per-sample directories with Cell Ranger-compatible naming
    (<sample_name>_S1_L00N_{R1|R2|I1|I2}_001.fastq.gz). sample_name is the ENA sample_title
    (sanitised). Every lane must have exactly one R1 and one R2, otherwise the rule fails and
    groups nothing.
    Output: {out_location}/ena_download/{project}/grouped/{sample_name}/...
    '''
    input:
        metadata  = rules.ena_fetch_metadata.output.metadata,
        downloads = get_ena_downloads
    output:
        done = touch(ENA_OUT + "{project}/grouped/.done")
    params:
        fastq_dir = ENA_OUT + "{project}/fastq",
        out_dir   = ENA_OUT + "{project}/grouped"
    wildcard_constraints:
        project = r"[A-Za-z0-9]+"
    conda:
        ENA["env_metadata"]
    log:
        "logs/ena_group_fastqs/{project}.log"
    benchmark:
        "benchmarks/ena_group_fastqs/{project}.tsv"
    script:
        "../scripts/group_ena_fastqs.py"


# ---------------------------------------------------------------------------- #
rule ena_make_samplesheet:
    '''
    Generate a samplesheet CSV from the grouped FASTQ directories (same script and format
    as make_samplesheet of the geo_download module).
    Columns: sample_name, sample_run, sample_path.
    '''
    input:
        done = rules.ena_group_fastqs.output.done
    output:
        csv = ENA_OUT + "{project}/samplesheet.csv"
    params:
        grouped_dir = ENA_OUT + "{project}/grouped"
    wildcard_constraints:
        project = r"[A-Za-z0-9]+"
    conda:
        ENA["env_metadata"]
    log:
        "logs/ena_make_samplesheet/{project}.log"
    benchmark:
        "benchmarks/ena_make_samplesheet/{project}.tsv"
    script:
        "../scripts/make_samplesheet.py"
