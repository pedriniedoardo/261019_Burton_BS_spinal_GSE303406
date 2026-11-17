"""
Fetch SRR run accessions and metadata for a GSE, BioProject (PRJNA/PRJEB/PRJDB) or
SRA study (SRP/ERP/DRP) accession using pysradb.
Called as a Snakemake script; accesses snakemake.params, .output, .log.
"""

import re
import sys
import logging

logging.basicConfig(
    level=logging.INFO,
    filename=snakemake.log[0],
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

try:
    import pandas as pd
    from pysradb.sraweb import SRAweb

    gse      = snakemake.params.gse
    out_srr  = snakemake.output.srr_list
    out_meta = snakemake.output.metadata

    logging.info(f"Fetching metadata for {gse}")

    db = SRAweb()

    if re.fullmatch(r"GSE\d+", gse):
        # GSE → SRP (a single GSE can map to multiple SRPs in rare cases)
        srp_df = db.gse_to_srp(gse)
        if srp_df is None or srp_df.empty:
            raise ValueError(f"No SRP found for {gse}. Check that the accession is correct.")
        srps = srp_df["study_accession"].unique().tolist()
        logging.info(f"Found SRP(s): {srps}")
    elif re.fullmatch(r"PRJ(NA|EB|DB)\d+|[SED]RP\d+", gse):
        # BioProject / SRA study: no GEO mapping needed, query SRA directly
        srps = [gse]
    else:
        raise ValueError(
            f"Unrecognised accession '{gse}': expected GSE<n>, PRJNA/PRJEB/PRJDB<n> or SRP/ERP/DRP<n>"
        )

    # SRP(s) → SRR
    all_runs = []
    for srp in srps:
        try:
            run_df = db.srp_to_srr(srp, detailed=True)
        except Exception as exc:
            logging.warning(f"srp_to_srr failed for {srp} ({exc}); trying sra_metadata")
            try:
                run_df = db.sra_metadata(srp, detailed=True)
            except Exception as exc2:
                # pysradb fails (e.g. AttributeError) when SRA has no runs for the study
                logging.warning(f"sra_metadata failed for {srp} ({exc2})")
                run_df = None
        if run_df is not None and not run_df.empty:
            all_runs.append(run_df)
        else:
            logging.warning(f"No SRR found for {srp}, skipping")

    if not all_runs:
        raise ValueError(
            f"SRA has no public runs for {gse} ({', '.join(srps)}). The reads may be hosted elsewhere: "
            f"check https://www.ebi.ac.uk/ena/browser/view/{srps[0]} and, if they are in ENA, "
            f"use `snakemake ena_samplesheet` (config['ena_download']['project'])."
        )

    combined = pd.concat(all_runs).drop_duplicates(subset="run_accession").reset_index(drop=True)
    srr_list = combined["run_accession"].tolist()
    logging.info(f"Found {len(srr_list)} SRR accession(s): {srr_list}")

    # columns used downstream for sample naming; absent for some non-GEO projects
    for col in ("experiment_title", "library_name", "sample_alias"):
        if col not in combined.columns:
            combined[col] = ""

    combined.to_csv(out_meta, index=False, sep="\t")

    with open(out_srr, "w") as fh:
        fh.write("\n".join(srr_list) + "\n")

    logging.info("Done")

except Exception as exc:
    logging.error(str(exc))
    sys.exit(1)
