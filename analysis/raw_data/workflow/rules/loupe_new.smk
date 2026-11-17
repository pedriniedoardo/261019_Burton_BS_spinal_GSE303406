# Rules to produce cloupe and h5ad objects from rds R Seurat objects

import os

# ── helper: estimate memory from RDS file size ────────────────────────────────
# scRNA-seq RDS files are typically gzip-compressed; in-memory expansion
# is ~10–20×. We use 15× + a 6 GB fixed overhead, floored at 32 GB.
def mem_from_rds(wildcards, input):
    rds_mb = os.path.getsize(input.rds) / 1e6      # compressed MB
    return max(32_000, int(rds_mb * 15) + 6_000)   # result in MB


rule seurat_to_cloupe_h5ad:
    input:
        rds = rules.addAnnotations.output.rds
    output:
        cloupe = config["out_location"] + config["out_location_seurat"] + "object/integrated_obj_" + CELLBENDER_TAG + "_annotated.cloupe",
        h5ad   = config["out_location"] + config["out_location_seurat"] + "object/integrated_obj_" + CELLBENDER_TAG + "_annotated.h5ad"
    log:
        "logs/seurat_to_cloupe_h5ad/rds2loupe_h5ad.log"
    benchmark:
        "benchmarks/Seurat/integrated/rds2loupe_h5ad.txt"
    params:
        umap_name = "umap",   # name of the Harmony-UMAP reduction in the object
    threads: 2
    resources:
        mem_mb = mem_from_rds,
    conda:
        "../envs/env_R_COSR_export_new.yaml"
    script:
        "../scripts/seurat_to_cloupe_h5ad_new.R"
