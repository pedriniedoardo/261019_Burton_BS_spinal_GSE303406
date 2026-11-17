"""
Group the downloaded ENA FASTQs per sample with Cell Ranger-compatible names.
Called as a Snakemake script; accesses snakemake.input, .params, .log.

Sample name = ENA sample_title (fallback: sample_alias, sample_accession) with every
character outside [A-Za-z0-9_.-] replaced by "_". The read type (I1/I2/R1/R2) and the
lane are taken from the submitted Cell Ranger-style file name
(<prefix>_S<n>_L<lane>_<type>_001.fastq.gz). Lanes are renumbered L001..L00N per sample.
Nothing is linked unless every lane has exactly one R1 and one R2.
"""

import csv
import re
import logging
from collections import defaultdict
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

FASTQ_RE = re.compile(r"_L(\d{3})_(I1|I2|R1|R2)_001\.fastq\.gz$")


def sanitise(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name.strip())


try:
    fastq_dir = Path(snakemake.params.fastq_dir)
    out_dir   = Path(snakemake.params.out_dir)

    with open(snakemake.input.metadata, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # run -> sample name
    run_to_name = {}
    for r in rows:
        raw = r["sample_title"] or r["sample_alias"] or r["sample_accession"]
        run_to_name[r["run_accession"]] = sanitise(raw)

    # A sample name must map to one ENA sample; several runs of the same sample are lanes
    name_to_samples = defaultdict(set)
    for r in rows:
        name_to_samples[run_to_name[r["run_accession"]]].add(r["sample_accession"])
    dup = {n: sorted(s) for n, s in name_to_samples.items() if len(s) > 1}
    if dup:
        raise ValueError(f"Sample names are not unique after sanitising: {dup}")

    # Pass 1: collect (sample, run, original lane) -> {read type: file}, check all of them
    lanes = defaultdict(dict)
    problems = []
    for r in rows:
        run, name = r["run_accession"], run_to_name[r["run_accession"]]
        for url in filter(None, r["submitted_ftp"].split(";")):
            src = fastq_dir / run / Path(url).name
            m = FASTQ_RE.search(src.name)
            if not m:
                problems.append(f"{run}: {src.name} is not named <prefix>_L00N_<I1|I2|R1|R2>_001.fastq.gz")
            elif not src.exists():
                problems.append(f"{run}: {src} is missing")
            else:
                lanes[(name, run, m.group(1))][m.group(2)] = src
    for (name, run, lane), files in lanes.items():
        for need in ("R1", "R2"):
            if need not in files:
                problems.append(f"{run} ({name}, L{lane}): no {need} file")
    if problems:
        raise ValueError("Cannot group FASTQs, nothing linked:\n  " + "\n  ".join(problems))

    # Pass 2: renumber lanes per sample and write the links
    per_sample = defaultdict(list)
    for key in sorted(lanes):
        per_sample[key[0]].append(key)

    logging.info(f"Found {len(per_sample)} sample(s)")
    for name, keys in sorted(per_sample.items()):
        sample_dir = out_dir / name
        sample_dir.mkdir(parents=True, exist_ok=True)
        for lane_num, key in enumerate(keys, start=1):
            for read_type, src in sorted(lanes[key].items()):
                link = sample_dir / f"{name}_S1_L{lane_num:03d}_{read_type}_001.fastq.gz"
                if link.is_symlink() or link.exists():
                    link.unlink()
                link.symlink_to(src.resolve())
                logging.info(f"  {src.name} -> {name}/{link.name}")

    logging.info("Done")

except Exception as exc:
    logging.exception(str(exc))
    raise
