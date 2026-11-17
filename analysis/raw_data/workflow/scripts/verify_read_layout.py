"""
Decide the final R1/R2/I1/I2 label of each FASTQ of one SRR and save it as JSON:
  [{"index": 0, "file": "SRRxxx_1.fastq.gz", "label": "R1", "source": "...",
    "read_length": 26, "barcode_match": 0.94, "whitelist": "737K-august-2016 (3' v2 / 5' v1-v2)"}, ...]

The labels are checked against the reads themselves: the first 16 bp of the first `n_reads` reads
of each file are compared with the 10x barcode whitelists shipped with the cellranger the pipeline
runs (locate_barcode_whitelists). The barcode read (R1) is the one file whose reads mostly start with
a valid barcode. Decision order:
  1. config override (read_layout_override; key SRR > GSM > GSE)
       - length must match the number of FASTQs
       - SRR/GSM key: used even if the barcodes disagree (WARNING)
       - GSE key:     ERROR if the barcodes point at another R1 (samples of a GSE may differ)
  2. SRA metadata layout (get_sra_read_layout), if complete (every file labelled, one R1 + one R2)
       - barcodes agree          → '<source>+verified'
       - barcodes disagree       → ERROR
       - no barcode read found   → '<source>+unverified' (WARNING; e.g. non-10x data)
  3. barcodes only: barcode read → R1, reads < 20 bp → I1, I2 (file order), the one remaining read → R2
  4. read length, only if allow_length_fallback and the pattern is unambiguous (WARNING)
  otherwise ERROR: set geo_download.read_layout_override for this SRR/GSM/GSE.
"""

import csv
import gzip
import json
import logging
import re
import statistics
from collections import Counter
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

BARCODE_LEN = 16
INDEX_MAX_LEN = 20    # reads shorter than this are sample-index reads (I1, I2)
CDNA_MIN_LEN = 50     # length fallback: >= this is cDNA, [INDEX_MAX_LEN, CDNA_MIN_LEN) is barcode+UMI

# Whitelist file stems in <cellranger>/lib/python/cellranger/barcodes (names differ slightly across
# cellranger versions, e.g. 3M-february-2018.txt.gz in 8.x vs 3M-february-2018_TRU.txt.gz in 10.x)
WHITELISTS = [
    ("737K-april-2014_rc",       "3' v1"),
    ("737K-august-2016",         "3' v2 / 5' v1-v2"),
    ("3M-february-2018",         "3' v3"),
    ("3M-3pgex-may-2023",        "3' v4"),
    ("3M-5pgex-jan-2023",        "5' v3"),
    ("737K-arc-v1",              "Multiome GEX"),
    ("737K-fixed-rna-profiling", "Flex"),
    ("737K-flex-v2",             "Flex v2"),
    ("9K-LT-march-2021",         "3' LT"),
    ("visium-v1",                "Visium v1"),
    ("visium-v2",                "Visium v2"),
    ("visium-v3",                "Visium v3"),
    ("visium-v4",                "Visium v4"),
    ("visium-v5",                "Visium v5"),
]

# Depositor labels that differ from the standard R1/R2/I1/I2 convention (from group_fastqs.py)
_LABEL_MAP = {
    "forward": "R1", "1": "R1", "read1": "R1", "r1": "R1",
    "reverse": "R2", "2": "R2", "read2": "R2", "r2": "R2",
    "index1":  "I1", "i1": "I1",
    "index2":  "I2", "i2": "I2",
    "barcode": "R1",
}


def normalise(label):
    raw = str(label or "").strip()
    return _LABEL_MAP.get(raw.lower(), raw.upper())


def list_fastqs(srr_dir, srr):
    """Return {0-based SRA read index: path}; fasterq-dump --split-files writes SRR_1, SRR_2, ..."""
    files = {}
    for f in srr_dir.glob(f"{srr}_*.fastq.gz"):
        m = re.search(r"_(\d+)\.fastq\.gz$", f.name)
        if m:
            files[int(m.group(1)) - 1] = f
    if not files and (srr_dir / f"{srr}.fastq.gz").exists():
        files[0] = srr_dir / f"{srr}.fastq.gz"
    return dict(sorted(files.items()))


def sample_reads(path, n_reads):
    """Return (median read length, Counter of first-16-bp prefixes, reads sampled)."""
    lengths, prefixes = [], Counter()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seq = line.strip()
                lengths.append(len(seq))
                if len(seq) >= BARCODE_LEN:
                    prefixes[seq[:BARCODE_LEN]] += 1
                if len(lengths) >= n_reads:
                    break
    if not lengths:
        return 0, prefixes, 0
    return int(statistics.median(lengths)), prefixes, len(lengths)


def find_whitelists(barcode_dir):
    found = []
    for stem, chemistry in WHITELISTS:
        for suffix in (".txt", ".txt.gz", "_TRU.txt", "_TRU.txt.gz"):
            p = barcode_dir / f"{stem}{suffix}"
            if p.exists():
                found.append((stem, chemistry, p))
                break
    return found


def whitelist_hits(path, wanted):
    """Stream one whitelist and return the subset of `wanted` prefixes it contains."""
    hits = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            bc = line[:BARCODE_LEN]
            if bc in wanted:
                hits.add(bc)
    return hits


def barcode_scan(reads, whitelists):
    """For each file: best (rate, whitelist name) over all whitelists."""
    wanted = set()
    for r in reads.values():
        wanted.update(r["prefixes"])
    best = {idx: (0.0, None) for idx in reads}
    for stem, chemistry, path in whitelists:
        hits = whitelist_hits(path, wanted)
        for idx, r in reads.items():
            if not r["n"]:
                continue
            rate = sum(c for bc, c in r["prefixes"].items() if bc in hits) / r["n"]
            if rate > best[idx][0]:
                best[idx] = (rate, f"{stem} ({chemistry})")
    return best


def gsm_for_srr(metadata_tsv, srr):
    with open(metadata_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("run_accession") == srr:
                alias = (row.get("experiment_alias") or "").strip()
                if alias.startswith("GSM"):
                    return alias
                m = re.match(r"(GSM\d+)", row.get("experiment_title") or "")
                return m.group(1) if m else None
    return None


def labels_by_barcode(reads, barcode_idx):
    """Barcode read → R1, short reads → I1, I2 (file order), the single remaining read → R2."""
    labels = {barcode_idx: "R1"}
    short = [i for i in reads if i != barcode_idx and reads[i]["length"] < INDEX_MAX_LEN]
    long_ = [i for i in reads if i != barcode_idx and reads[i]["length"] >= INDEX_MAX_LEN]
    if len(long_) != 1:
        return None
    labels[long_[0]] = "R2"
    for n, i in enumerate(short, start=1):
        labels[i] = f"I{n}"
    return labels


def labels_by_length(reads):
    """Only the unambiguous 10x pattern: one 20-49 bp read, one >= 50 bp read, the rest < 20 bp."""
    mid = [i for i, r in reads.items() if INDEX_MAX_LEN <= r["length"] < CDNA_MIN_LEN]
    long_ = [i for i, r in reads.items() if r["length"] >= CDNA_MIN_LEN]
    short = [i for i, r in reads.items() if r["length"] < INDEX_MAX_LEN]
    if len(mid) != 1 or len(long_) != 1:
        return None
    labels = {mid[0]: "R1", long_[0]: "R2"}
    for n, i in enumerate(short, start=1):
        labels[i] = f"I{n}"
    return labels


def complete(labels, reads):
    values = list(labels.values())
    return (set(labels) == set(reads) and len(values) == len(set(values))
            and values.count("R1") == 1 and values.count("R2") == 1)


def fmt(labels, reads):
    return ", ".join(f"{reads[i]['file'].name}={labels.get(i, '?')}" for i in reads)


try:
    srr       = snakemake.wildcards.srr
    gse       = snakemake.wildcards.gse
    srr_dir   = Path(snakemake.params.fastq_dir) / srr
    n_reads   = int(snakemake.params.n_reads)
    min_match = float(snakemake.params.min_barcode_match)
    allow_len = bool(snakemake.params.allow_length_fallback)
    overrides = {str(k): v for k, v in (snakemake.params.overrides or {}).items()}

    barcode_dir = Path(Path(snakemake.input.whitelist_dir).read_text().strip())
    whitelists = find_whitelists(barcode_dir)
    if not whitelists:
        raise RuntimeError(f"No 10x barcode whitelists found in {barcode_dir}")
    logging.info(f"{srr}: whitelists from {barcode_dir} ({len(whitelists)} lists)")

    fastqs = list_fastqs(srr_dir, srr)
    if not fastqs:
        raise RuntimeError(f"No FASTQ files found for {srr} in {srr_dir}")

    # 1. sample the reads and match them against the whitelists
    reads = {}
    for idx, f in fastqs.items():
        length, prefixes, n = sample_reads(f, n_reads)
        reads[idx] = {"file": f, "length": length, "prefixes": prefixes, "n": n}
    best = barcode_scan(reads, whitelists)
    logging.info(f"Barcode check on the first {n_reads} reads (threshold {min_match:.0%}):")
    for idx, r in reads.items():
        rate, wl = best[idx]
        logging.info(f"  {r['file'].name}: {r['length']} bp, {rate:.1%} valid barcodes"
                     + (f" (best: {wl})" if wl else ""))

    matching = [i for i in reads if best[i][0] >= min_match]
    barcode_idx = matching[0] if len(matching) == 1 else None
    if len(matching) > 1:
        logging.warning(f"  more than one file matches the whitelists: "
                        f"{[reads[i]['file'].name for i in matching]}; barcode read not identified")
    elif barcode_idx is not None:
        logging.info(f"  barcode read (R1): {reads[barcode_idx]['file'].name}, {best[barcode_idx][1]}")
    else:
        logging.info("  no file matches a 10x barcode whitelist")

    def barcode_disagrees(labels):
        r1 = [i for i, l in labels.items() if l == "R1"]
        return barcode_idx is not None and r1 != [barcode_idx]

    # 2. decide
    labels, source = None, None
    gsm = gsm_for_srr(snakemake.input.metadata, srr)
    level, override = next(((lvl, overrides[key]) for lvl, key in (("srr", srr), ("gsm", gsm), ("gse", gse))
                            if key and key in overrides), (None, None))
    keys = " or ".join(k for k in (srr, gsm) if k)   # per-sample override keys for error messages
    meta_layout = json.loads(Path(snakemake.input.metadata_layout).read_text())
    meta = {r["index"]: normalise(r["label"]) for r in meta_layout if r.get("label")}
    meta_source = meta_layout[0].get("source", "metadata") if meta_layout else "metadata"

    if override is not None:
        override = [normalise(l) for l in override]
        key = {"srr": srr, "gsm": gsm, "gse": gse}[level]
        if len(override) != len(reads):
            raise RuntimeError(
                f"read_layout_override[{key}] has {len(override)} labels but {srr} has "
                f"{len(reads)} FASTQs ({', '.join(r['file'].name for r in reads.values())})")
        labels = dict(zip(reads, override))
        source = f"config-override-{level}"
        if barcode_disagrees(labels):
            msg = (f"read_layout_override[{key}] ({fmt(labels, reads)}) disagrees with the barcode check, "
                   f"which found R1 = {reads[barcode_idx]['file'].name}")
            if level == "gse":
                raise RuntimeError(msg + f". Samples of a GSE may differ: add a per-sample entry ({keys}) to read_layout_override")
            logging.warning(msg + f"; using the {level.upper()} override as requested")
        elif barcode_idx is None:
            logging.warning("  override not verified: no barcode read found")

    elif meta and complete(meta, reads):
        labels = meta
        if barcode_disagrees(labels):
            raise RuntimeError(
                f"SRA metadata ({meta_source}: {fmt(labels, reads)}) disagrees with the barcode check, "
                f"which found R1 = {reads[barcode_idx]['file'].name} ({best[barcode_idx][0]:.1%}). "
                f"Check the reads and set read_layout_override for {keys}")
        if barcode_idx is None:
            source = f"{meta_source}+unverified"
            logging.warning("  metadata labels kept but not verified (no 10x barcode read found)")
        else:
            source = f"{meta_source}+verified"

    elif barcode_idx is not None:
        if meta:
            logging.warning(f"  metadata layout incomplete ({fmt(meta, reads)}); using the barcode check")
        labels = labels_by_barcode(reads, barcode_idx)
        if labels is None:
            raise RuntimeError(
                f"Barcode read is {reads[barcode_idx]['file'].name}, but the other reads do not contain "
                f"exactly one cDNA read (>= {INDEX_MAX_LEN} bp). Set read_layout_override for {keys}")
        source = "barcode-whitelist"

    elif allow_len and (labels := labels_by_length(reads)) is not None:
        source = "length-heuristic"
        logging.warning(f"  labels guessed from read length (allow_length_fallback): {fmt(labels, reads)}")

    else:
        raise RuntimeError(
            f"Could not assign R1/R2 for {srr}" + (f" ({gsm})" if gsm else "") + ": no usable SRA layout, no read matches a 10x "
            f"barcode whitelist" + ("" if allow_len else ", length fallback disabled")
            + f". Check the reads and set geo_download.read_layout_override for {keys} or {gse}")

    layout = [
        {
            "index": idx,
            "file": reads[idx]["file"].name,
            "label": labels[idx],
            "source": source,
            "read_length": reads[idx]["length"],
            "barcode_match": round(best[idx][0], 4),
            "whitelist": best[idx][1] if best[idx][0] >= min_match else None,
        }
        for idx in reads
    ]
    logging.info(f"Final layout ({source}): {fmt(labels, reads)}")
    with open(snakemake.output.layout, "w") as fh:
        json.dump(layout, fh, indent=2)

except Exception as exc:
    logging.error(str(exc))
    raise
