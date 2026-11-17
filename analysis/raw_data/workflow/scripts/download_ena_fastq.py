"""
Download the submitted FASTQ files of one ENA run and verify their MD5.
Called as a Snakemake script; accesses snakemake.input, .params, .log.

Resumable: data goes to <file>.part and is resumed with an HTTP Range request;
a file that already exists with the right MD5 is not downloaded again.
A connection dropped mid-file is detected by comparing the received bytes with the size the
server announced (http.client does not raise on a short read(amt)); the .part is kept and
resumed. ENA throttles clients (HTTP 403/429/503, refused connections), so those errors are
retried with a long, growing pause. Only a complete-size file with a wrong MD5 is discarded.
"""

import csv
import re
import sys
import time
import hashlib
import logging
import urllib.error
import urllib.request
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

CHUNK = 1 << 20


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


THROTTLE_CODES = {403, 429, 503}


class Truncated(Exception):
    pass


def expected_total(resp, start):
    """Total file size announced by the server (None if unknown)."""
    cr = resp.headers.get("Content-Range")            # "bytes 100-999/1000"
    if cr:
        m = re.match(r"bytes \d+-\d+/(\d+)$", cr)
        return int(m.group(1)) if m else None
    cl = resp.headers.get("Content-Length")
    return start + int(cl) if cl is not None else None


def fetch(url, dest, retries=8, max_attempts=500):
    """Download url to dest.part (resuming), return the .part path once it has the announced size.

    Only attempts that add no bytes count towards `retries`: ENA often closes the connection
    every ~1 GB, and such an attempt still moves the download forward. `max_attempts` bounds
    the total number of attempts.
    """
    part = Path(str(dest) + ".part")
    stalled = 0                                     # consecutive attempts without progress
    for attempt in range(1, max_attempts + 1):
        start = part.stat().st_size if part.exists() else 0
        req = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                if start and resp.status != 206:   # server ignored the Range header: restart
                    start = 0
                total = expected_total(resp, start)
                with open(part, "ab" if start else "wb") as fh:
                    for block in iter(lambda: resp.read(CHUNK), b""):
                        fh.write(block)
            size = part.stat().st_size
            if total is not None and size != total:
                raise Truncated(f"connection closed at {size} of {total} bytes")
            return part
        except urllib.error.HTTPError as exc:
            if exc.code == 416:                     # range not satisfiable: .part is at/after the end
                return part
            reason, throttled = f"HTTP {exc.code}", exc.code in THROTTLE_CODES
        except Truncated as exc:
            reason, throttled = str(exc), False
        except Exception as exc:
            reason, throttled = str(exc), True      # refused/reset/timeout: usually throttling too
        size = part.stat().st_size if part.exists() else 0
        stalled = 0 if size > start else stalled + 1
        if stalled >= retries:
            raise RuntimeError(f"Download failed, no progress in {retries} attempts: {url} ({reason})")
        pause = (30 if throttled else 5) * max(stalled, 1)
        logging.warning(f"{url}: {reason}; at {size} bytes, retrying in {pause}s "
                        f"(attempt {attempt}, {stalled}/{retries} without progress)")
        time.sleep(pause)
    raise RuntimeError(f"Download failed after {max_attempts} attempts: {url}")


try:
    run    = snakemake.params.run
    outdir = Path(snakemake.params.outdir)

    with open(snakemake.input.metadata, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["run_accession"] == run]
    if len(rows) != 1:
        raise ValueError(f"Expected 1 metadata row for {run}, found {len(rows)}")
    row = rows[0]

    urls = [u for u in row["submitted_ftp"].split(";") if u]
    md5s = [m for m in row["submitted_md5"].split(";") if m]
    outdir.mkdir(parents=True, exist_ok=True)

    for url, expected in zip(urls, md5s):
        dest = outdir / Path(url).name
        if dest.exists() and md5sum(dest) == expected:
            logging.info(f"{dest.name}: already downloaded, md5 OK")
            continue
        https = url if url.startswith(("http://", "https://")) else "https://" + url
        logging.info(f"{dest.name}: downloading {https}")
        part = fetch(https, dest)
        got = md5sum(part)
        if got != expected:
            part.unlink()   # full-size but corrupt: do not resume from it next time
            raise ValueError(f"{dest.name}: md5 mismatch (expected {expected}, got {got})")
        part.rename(dest)
        logging.info(f"{dest.name}: md5 OK")

    logging.info(f"Done: {len(urls)} file(s) for {run}")

except Exception as exc:
    logging.exception(str(exc))
    sys.exit(1)
