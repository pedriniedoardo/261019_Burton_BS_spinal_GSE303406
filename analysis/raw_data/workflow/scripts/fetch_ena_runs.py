"""
Fetch the run-level metadata of an ENA/SRA project from the ENA portal API.
Called as a Snakemake script; accesses snakemake.params, .output, .log.

Only projects whose submitted files are FASTQs are supported: the FASTQs are
downloaded as submitted (see download_ena_fastq.py / group_ena_fastqs.py).
"""

import csv
import io
import sys
import time
import logging
import urllib.parse
import urllib.request

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

API = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = [
    "run_accession", "experiment_accession", "sample_accession", "sample_alias",
    "sample_title", "library_layout", "submitted_format", "submitted_ftp", "submitted_md5",
]


def query_ena(project, retries=3):
    params = {
        "result": "read_run",
        "query": f'(study_accession="{project}" OR secondary_study_accession="{project}")',
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": 0,
    }
    url = API + "?" + urllib.parse.urlencode(params)
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                return resp.read().decode("utf-8")
        except Exception as exc:
            logging.warning(f"ENA query failed (attempt {attempt}/{retries}): {exc}")
            if attempt == retries:
                raise
            time.sleep(5 * attempt)


try:
    project  = snakemake.params.project
    out_runs = snakemake.output.run_list
    out_meta = snakemake.output.metadata

    logging.info(f"Fetching ENA runs for {project}")
    text = query_ena(project)

    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise ValueError(
            f"No public runs found in ENA for {project}. Check the accession at "
            f"https://www.ebi.ac.uk/ena/browser/view/{project}"
        )

    problems = []
    for r in rows:
        run = r["run_accession"]
        files = [x for x in r["submitted_ftp"].split(";") if x]
        md5s  = [x for x in r["submitted_md5"].split(";") if x]
        formats = {x for x in r["submitted_format"].split(";") if x}
        if not files:
            problems.append(f"{run}: no submitted files (only ENA-generated FASTQs?)")
        elif formats != {"FASTQ"}:
            problems.append(f"{run}: submitted files are {sorted(formats)}, only FASTQ is supported")
        elif len(md5s) != len(files):
            problems.append(f"{run}: {len(files)} files but {len(md5s)} md5 checksums")
    if problems:
        raise ValueError("Unsupported ENA project:\n  " + "\n  ".join(problems))

    runs = sorted({r["run_accession"] for r in rows})
    logging.info(f"Found {len(runs)} run(s), "
                 f"{len({r['sample_accession'] for r in rows})} sample(s): {runs}")

    with open(out_meta, "w") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    with open(out_runs, "w") as fh:
        fh.write("\n".join(runs) + "\n")

    logging.info("Done")

except Exception as exc:
    logging.exception(str(exc))
    sys.exit(1)
