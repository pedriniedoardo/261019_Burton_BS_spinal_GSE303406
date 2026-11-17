"""
Build a per-SRR read-label mapping and save it as JSON:
  [{"index": 0, "label": "I1", "source": "efetch"},  ...]

Three-tier approach:
  1. efetch READ_SPEC  — parses READ_LABEL from the SRA SPOT_DESCRIPTOR (most authoritative).
                         Present for well-annotated submissions.
  2. efetch SRAFile    — parses original submitted filenames from SRAFiles section.
                         Filenames like *_R1_001.fastq.gz, *_I1_001.fastq.gz encode read type
                         directly. Works for GEO-submitted data that skipped SPOT_DESCRIPTOR.
  3. sra-stat          — fallback: infers labels from Technical/Biological + read length.
  4. empty list        — if all fail.
This is metadata only: verify_read_layout.py checks it against the reads (10x barcode whitelists)
and writes the final layout used by group_fastqs.py.
"""

import json
import logging
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


_READ_TYPE_RE = re.compile(r'[_\-]([IR][12])[_\-]001\.fastq', re.IGNORECASE)

def _fetch_xml(srr):
    """Run efetch and return the parsed XML tree, or None on failure."""
    try:
        result = subprocess.run(
            ["efetch", "-db", "sra", "-id", srr, "-format", "xml"],
            capture_output=True, text=True, check=True, timeout=120
        )
        return ET.fromstring(result.stdout)
    except Exception as e:
        logging.warning(f"efetch failed: {e}")
        return None


def try_efetch_read_spec(tree):
    """Tier 1: parse READ_LABEL from SPOT_DESCRIPTOR/READ_SPEC (most authoritative)."""
    if tree is None:
        return []
    specs = [
        {
            "index":  int(spec.findtext("READ_INDEX")),
            "label":  spec.findtext("READ_LABEL") or "",
            "source": "efetch-READ_SPEC"
        }
        for spec in tree.findall(".//READ_SPEC")
        if spec.findtext("READ_LABEL")
    ]
    return specs


def try_efetch_filenames(tree, srr):
    """
    Tier 2: parse read type from original submitted filenames in SRAFiles section.
    Files with supertype='Original' have names like *_R1_001.fastq.gz, *_I1_001.fastq.gz.
    Their order in the XML matches fasterq-dump's _1, _2, _3 output order (READ_INDEX 0-based).
    Only the RUN matching `srr` is read: efetch returns the whole experiment package, which can
    list sibling runs of the same sample, whose files would otherwise shift the indices.
    """
    if tree is None:
        return []
    runs = [r for r in tree.iter("RUN") if r.get("accession") == srr]
    if not runs:
        logging.warning(f"  {srr} not found among the RUNs returned by efetch")
        return []
    originals = [
        f.get("filename", "")
        for f in runs[0].findall(".//SRAFile")
        if f.get("supertype") == "Original" and f.get("filename", "").endswith(".fastq.gz")
    ]
    if not originals:
        return []
    result = []
    for idx, fname in enumerate(originals):
        m = _READ_TYPE_RE.search(fname)
        if m:
            result.append({"index": idx, "label": m.group(1).upper(), "source": "efetch-filename"})
        else:
            logging.warning(f"  Could not extract read type from original filename: {fname}")
    return result


def try_sra_stat(sra_path):
    """Return list of {index, label, source} inferred from sra-stat Technical/Biological."""
    try:
        result = subprocess.run(
            ["sra-stat", "--quick", "--xml", sra_path],
            capture_output=True, text=True, check=True, timeout=120
        )
        tree = ET.fromstring(result.stdout)
        reads = [
            {
                "index":      int(el.get("index")),
                "type":       el.get("type"),           # "Technical" or "Biological"
                "avg_length": float(el.get("average"))
            }
            for el in tree.findall(".//Read")
        ]
        if not reads:
            return []

        biological = sorted([r for r in reads if r["type"] == "Biological"], key=lambda r: r["index"])
        technical  = sorted([r for r in reads if r["type"] == "Technical"],  key=lambda r: r["avg_length"])

        result_list = []

        # Biological → R2 (cDNA); edge case: no technical reads → single-end R1
        if biological and not technical:
            result_list.append({"index": biological[0]["index"], "label": "R1", "source": "sra-stat"})
        else:
            for i, r in enumerate(biological, start=2):
                result_list.append({"index": r["index"], "label": f"R{i}", "source": "sra-stat"})

        # Technical: longest → R1 (barcode+UMI), shorter → I1, I2
        index_labels = ["I1", "I2", "I3", "I4"]
        for i, r in enumerate(technical[:-1]):
            result_list.append({"index": r["index"], "label": index_labels[i], "source": "sra-stat"})
        if technical:
            result_list.append({"index": technical[-1]["index"], "label": "R1", "source": "sra-stat"})

        return sorted(result_list, key=lambda r: r["index"])

    except Exception as e:
        logging.warning(f"sra-stat failed: {e}")
        return []


try:
    srr      = Path(snakemake.input.sra).parent.name
    sra_path = snakemake.input.sra
    out      = snakemake.output.layout

    xml_tree = _fetch_xml(srr)

    layout = try_efetch_read_spec(xml_tree)
    if layout:
        logging.info(f"Source: efetch READ_SPEC ({len(layout)} reads)")
    else:
        logging.info("No READ_SPEC — trying original filenames from SRAFiles")
        layout = try_efetch_filenames(xml_tree, srr)
        if layout:
            logging.info(f"Source: efetch original filenames ({len(layout)} reads)")
        else:
            logging.info("No original filenames — trying sra-stat")
            layout = try_sra_stat(sra_path)
            if layout:
                logging.info(f"Source: sra-stat inferred ({len(layout)} reads)")
            else:
                logging.warning("All sources failed — writing empty layout; verify_read_layout.py will detect it from the reads")

    for r in layout:
        logging.info(f"  index={r['index']}  label='{r['label']}'  source={r['source']}")

    with open(out, "w") as fh:
        json.dump(layout, fh, indent=2)

except Exception as exc:
    logging.error(str(exc))
    raise
