"""
Generate a samplesheet CSV from the per-sample grouped FASTQ directories.
Output columns: sample_name, sample_run, sample_path
sample_path is the absolute path to the directory — compatible with the
SpaceRanger pipeline's samples_csv input format.
"""

import logging
import csv
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

try:
    grouped_dir = Path(snakemake.params.grouped_dir).resolve()
    out_csv     = snakemake.output.csv

    sample_dirs = sorted(
        d for d in grouped_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    )

    if not sample_dirs:
        raise ValueError(f"No sample directories found in {grouped_dir}")

    rows = [
        {
            "sample_name": d.name,
            "sample_run":  "run01",
            "sample_path": str(d)
        }
        for d in sample_dirs
    ]

    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_name", "sample_run", "sample_path"])
        writer.writeheader()
        writer.writerows(rows)

    logging.info(f"Written {len(rows)} sample(s) to {out_csv}")
    for r in rows:
        logging.info(f"  {r['sample_name']}  {r['sample_run']}  {r['sample_path']}")

except Exception as exc:
    logging.error(str(exc))
    raise
