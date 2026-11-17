# ======================================================================
# == load libraries ==
# ======================================================================

# loupeR: install if absent
if (!requireNamespace("loupeR", quietly = TRUE)) {
  remotes::install_github("10XGenomics/loupeR")
}

# Keep the louper binary inside the active conda env so all cluster nodes share it
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix == "") stop("CONDA_PREFIX not set — is this running inside a conda env?")
target_dir <- file.path(conda_prefix, "share", "loupeR")
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
Sys.setenv(LOUPER_USER_DATA_DIR = target_dir)

# Non-interactive EULA acceptance confirmed with 10X Genomics support (2025-05)
Sys.setenv(AUTO_ACCEPT_EULA = "true")

suppressPackageStartupMessages({
  library(Seurat)
  library(loupeR)
  library(anndataR)
})

# Download the louper binary on first run; no-op if already present
loupeR::setup()

# ── Redirect stdout + stderr to snakemake log ────────────────────────────────
log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) })

ts <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")

# ── Inputs / outputs / params ────────────────────────────────────────────────
input_rds     <- snakemake@input[["rds"]]
output_cloupe <- snakemake@output[["cloupe"]]
output_h5ad   <- snakemake@output[["h5ad"]]
umap_name     <- if (!is.null(snakemake@params[["umap_name"]])) {
                   snakemake@params[["umap_name"]]
                 } else "umap"

# ── Load object ──────────────────────────────────────────────────────────────
ts("Loading RDS:", input_rds)
seu <- readRDS(input_rds)
ts(sprintf("Loaded: %d cells × %d features", ncol(seu), nrow(seu)))
ts("Assays   :", paste(names(seu@assays), collapse = ", "))
ts("Reductions:", paste(names(seu@reductions), collapse = ", "))
ts("Metadata cols:", paste(colnames(seu@meta.data), collapse = ", "))

# ── Ensure RNA is default assay ──────────────────────────────────────────────
if (!"RNA" %in% names(seu@assays)) stop("No 'RNA' assay found in the Seurat object.")
DefaultAssay(seu) <- "RNA"

if (!umap_name %in% names(seu@reductions)) {
  stop(sprintf(
    "Reduction '%s' not found. Available: %s",
    umap_name, paste(names(seu@reductions), collapse = ", ")
  ))
}

# =============================================================================
# 1.  .cloupe  ──  loupeR
# =============================================================================
ts("Building .cloupe with loupeR...")
cloupe_dir  <- dirname(output_cloupe)
cloupe_name <- tools::file_path_sans_ext(basename(output_cloupe))

loupeR::create_loupe_from_seurat(
  obj         = seu,
  output_dir  = cloupe_dir,
  output_name = cloupe_name,
  force       = TRUE
)
ts("Cloupe written to:", output_cloupe)

# =============================================================================
# 2.  .h5ad  ──  anndataR (pure R, no Python / reticulate dependency)
# =============================================================================
ts("Writing .h5ad...")
write_h5ad(seu, output_h5ad)
ts("H5AD written to:", output_h5ad)
ts("All done.")
