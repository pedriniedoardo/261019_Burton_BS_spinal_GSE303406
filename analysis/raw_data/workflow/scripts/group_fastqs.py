"""
Group downloaded FASTQs by biological sample with SpaceRanger-compatible naming.
Each experiment (SRX) becomes one sample directory; multiple SRR runs become lanes (L001, L002, ...).

Read types come from read_layout/{srr}.json written by verify_read_layout.py (config override,
SRA metadata checked against the 10x barcode whitelists, barcode detection, or opt-in length
fallback). Before any link is made, every lane must have each FASTQ labelled, no duplicate labels
and exactly one R1 and one R2; otherwise nothing is grouped and all failing SRRs are reported.
"""

import re
import json
import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


# cellranger --id must be <= 64 characters and the sample name is reused for the FASTQ prefix,
# so keep names well below that
MAX_NAME_LEN = 50


def sanitize(name):
    return re.sub(r"[^\w\-]", "_", str(name)).strip("_")


def short_name(name):
    """
    Shorten a GEO experiment title such as
    'GSM4796271: LCL_777_B958; Homo sapiens; Human gammaherpesvirus 4; RNA-Seq'
    to 'GSM4796271_LCL_777_B958': keep the part before the first ';' (accession + sample label),
    sanitize it and collapse repeated underscores. The leading GSM/SRX accession keeps it unique.
    """
    first = str(name).split(";")[0]
    return re.sub(r"_+", "_", sanitize(first))[:MAX_NAME_LEN].strip("_")


def assign_read_types(fastq_files, read_layout):
    """
    Map each FASTQ to its label from the final layout written by verify_read_layout.py.
    fasterq-dump --split-files outputs _{index+1}.fastq.gz matching the layout's 0-based index.
    """
    def sra_index(f):
        m = re.search(r"_(\d+)\.fastq\.gz$", f.name)
        return int(m.group(1)) - 1 if m else 0   # filename 1-based → SRA 0-based

    labels = {r["index"]: r["label"] for r in read_layout}
    result = {}
    for f in fastq_files:
        label = labels.get(sra_index(f))
        result[f] = label
        logging.info(f"    {f.name}  label='{label}'")
    return result


def check_lane(srr, read_map):
    """Check A: every FASTQ labelled, no duplicate labels, exactly one R1 and one R2."""
    problems = []
    labels = list(read_map.values())
    unlabelled = [f.name for f, l in read_map.items() if not l]
    if unlabelled:
        problems.append(f"unlabelled {unlabelled}")
    dup = sorted({l for l in labels if l and labels.count(l) > 1})
    if dup:
        problems.append(f"duplicate labels {dup}")
    for req in ("R1", "R2"):
        if labels.count(req) != 1:
            problems.append(f"{labels.count(req)} x {req}")
    return [f"{srr}: " + "; ".join(problems)] if problems else []


try:
    meta      = pd.read_csv(snakemake.input.metadata, sep="\t")
    fastq_dir = Path(snakemake.params.fastq_dir)
    out_dir   = Path(snakemake.params.out_dir)

    # Build SRR → read layout mapping from the JSON files produced by get_sra_read_layout
    srr_to_layout = {}
    for lf in snakemake.input.layouts:
        srr = Path(lf).stem
        with open(lf) as fh:
            srr_to_layout[srr] = json.load(fh)

    if "experiment_accession" not in meta.columns:
        raise ValueError(
            "'experiment_accession' column missing from metadata.tsv — "
            "make sure fetch_srr_from_gse.py calls srp_to_srr(detailed=True)"
        )

    # SRR → SRX grouping (one SRX = one biological sample)
    groups = (
        meta.groupby("experiment_accession")["run_accession"]
        .apply(list)
        .to_dict()
    )

    # Pick best human-readable name (first non-empty column wins)
    name_priority = ["experiment_title", "library_name", "sample_alias", "experiment_accession"]
    srx_to_name = {}
    for srx, grp in meta.groupby("experiment_accession"):
        row = grp.iloc[0]
        for col in name_priority:
            if col in row.index and pd.notna(row[col]) and str(row[col]).strip():
                name = short_name(row[col])
                # titles of non-GEO projects carry no GSM accession: append the SRX to keep names unique
                if not re.match(r"GSM\d+", name) and srx not in name:
                    name = f"{name[:MAX_NAME_LEN - len(srx) - 1].strip('_')}_{srx}"
                srx_to_name[srx] = name
                break

    dup_names = pd.Series(srx_to_name).loc[lambda x: x.duplicated(keep=False)]
    if not dup_names.empty:
        raise ValueError(f"Sample names are not unique after shortening: {dup_names.to_dict()}")

    logging.info(f"Found {len(groups)} sample(s):")
    for srx, srrs in groups.items():
        logging.info(f"  {srx} → '{srx_to_name[srx]}': {srrs}")

    # Pass 1: assign labels to every lane and check them all before touching the output
    planned, problems = [], []
    for srx, srrs in groups.items():
        sample_name = srx_to_name[srx]
        logging.info(f"\nSample: {sample_name}")

        for lane_num, srr in enumerate(sorted(srrs), start=1):
            srr_dir = fastq_dir / srr
            fastqs  = sorted(srr_dir.glob(f"{srr}_*.fastq.gz")) or \
                      sorted(srr_dir.glob(f"{srr}.fastq.gz"))

            if not fastqs:
                problems.append(f"{srr}: no FASTQ files in {srr_dir}")
                continue

            layout = srr_to_layout.get(srr, [])
            sources = sorted({r.get("source", "?") for r in layout})
            logging.info(f"  {srr} → L{lane_num:03d}: {[f.name for f in fastqs]}  source={sources}")
            read_map = assign_read_types(fastqs, layout)
            problems += check_lane(srr, read_map)
            planned.append((sample_name, lane_num, read_map))

    if problems:
        raise ValueError(
            "Invalid read layout, nothing grouped:\n  " + "\n  ".join(problems)
            + "\nSee logs/verify_read_layout/ and set geo_download.read_layout_override if needed."
        )

    # Pass 2: all lanes valid → write the links
    for sample_name, lane_num, read_map in planned:
        sample_dir = out_dir / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        for src, read_type in read_map.items():
            link_name = f"{sample_name}_S1_L{lane_num:03d}_{read_type}_001.fastq.gz"
            link_path = sample_dir / link_name
            if link_path.is_symlink() or link_path.exists():
                link_path.unlink()
            link_path.symlink_to(src.resolve())
            logging.info(f"    {src.name} → {link_name}")

    logging.info("\nDone")

except Exception as exc:
    logging.error(str(exc))
    raise
