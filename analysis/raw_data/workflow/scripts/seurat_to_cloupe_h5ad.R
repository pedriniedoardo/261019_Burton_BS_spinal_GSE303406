# =============================================================================
# seurat_to_cloupe_h5ad.R
# Convert a Seurat v5 RDS object → .cloupe (Loupe Browser) + .h5ad (scanpy)
# Snakemake inputs:  rds
# Snakemake outputs: cloupe, h5ad
# Snakemake params:  umap_name  (default: "umap")
# =============================================================================


# install dependencies
if (!require("hdf5r")) install.packages("hdf5r")
if (!require("Seurat")) install.packages("Seurat")

library(remotes)
remotes::install_github("10XGenomics/loupeR")

suppressPackageStartupMessages({
  library(Seurat)
  library(loupeR)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Matrix)
})
# ── loupeR setup: keep everything inside the Snakemake conda env ─────────────
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix == "") stop("CONDA_PREFIX not set — is this running inside a conda env?")

# 1. Dynamically find the correct R user directory for the active user
target_dir <- file.path(conda_prefix, "share", "loupeR")
# 2. Force loupeR to strictly use this path (bypasses any cluster node environment wiping)
Sys.setenv(LOUPER_USER_DATA_DIR = target_dir)

# 3. If the directory doesn't exist yet for this user, create it
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)

# 4. Auto-approve the EULA file if it isn't there
eula_file <- file.path(target_dir, "eula_agreement")
if (!file.exists(eula_file)) writeLines("ACCEPTED", eula_file)

# 5. Load the library safely
library(loupeR)

# 6. If this user is missing the 32MB backend binary, download it headless

binary_file <- file.path(target_dir, "louper")
if (!file.exists(binary_file)) {
  message("louper binary missing — downloading into conda env...")
  
  louper_version <- as.character(packageVersion("loupeR"))
  os  <- Sys.info()[["sysname"]]
  bin <- if (os == "Darwin") "louper-macos-x64" else "louper-linux-x64"
  url <- paste0(
    "https://github.com/10XGenomics/loupeR/releases/download/v",
    louper_version, "/", bin
  )
  
  message("Downloading from: ", url)
  download.file(url, destfile = binary_file, mode = "wb", quiet = FALSE)
  Sys.chmod(binary_file, mode = "0755")
}

# NEVER call loupeR::setup() below this line# 1. Open a temporary text connection that holds the "yes" answer
# loupeR::setup()

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

# ── Seurat v5: join layers so counts are a single sparse matrix ───────────────
if (inherits(seu[["RNA"]], "Assay5")) {
  ts("Seurat v5 Assay5 detected — joining layers...")
  seu <- JoinLayers(seu)
  ts("Layers joined:", paste(Layers(seu[["RNA"]]), collapse = ", "))
}

# Sanity-check that the requested UMAP reduction exists
if (!umap_name %in% names(seu@reductions)) {
  stop(sprintf(
    "Reduction '%s' not found. Available: %s",
    umap_name, paste(names(seu@reductions), collapse = ", ")
  ))
}

# =============================================================================
# 1.  .cloupe  ──  loupeR
# =============================================================================
ts("Setting up loupeR (binary download on first run)...")
# loupeR::setup()

ts("Building .cloupe with loupeR...")
cloupe_dir  <- dirname(output_cloupe)
cloupe_name <- tools::file_path_sans_ext(basename(output_cloupe))
dir.create(cloupe_dir, recursive = TRUE, showWarnings = FALSE)

# create_cloupe_from_seurat uses:
#   • active assay counts  → RNA counts
#   • all @meta.data cols  → metadata / clusters
#   • all reductions       → projections (UMAP, PCA, …)
# We only need to make sure RNA is active and layers are joined (done above).
loupeR::create_loupe_from_seurat(
  obj         = seu,
  output_dir  = cloupe_dir,
  output_name = cloupe_name,
  force       = TRUE          # overwrite if rerunning
)
ts("Cloupe written to:", output_cloupe)

# =============================================================================
# 2.  .h5ad  ──  r-anndata (fully conda-native, no basilisk)
# =============================================================================
library(anndata)
library(reticulate)

# Use the conda env's own Python (activated by Snakemake)
use_python(Sys.which("python"), required = TRUE)

ts("Building AnnData object...")

# Cells × genes count matrix (AnnData convention: obs=cells, var=genes)
counts_mat <- t(GetAssayData(seu, assay = "RNA", layer = "counts"))

# Collect all reductions into obsm with X_ prefix (scanpy convention)
obsm_list <- list()
for (red in names(seu@reductions)) {
  key <- paste0("X_", gsub("[^A-Za-z0-9_]", "_", red))
  obsm_list[[key]] <- Embeddings(seu, reduction = red)
  ts(sprintf("  reduction '%s' → obsm '%s' (%d dims)", red, key, ncol(obsm_list[[key]])))
}

# var dataframe (one row per gene)
var_df <- data.frame(
  gene_name = rownames(seu),
  row.names  = rownames(seu)
)

adata <- AnnData(
  X    = counts_mat,
  obs  = seu@meta.data,
  var  = var_df,
  obsm = obsm_list
)

h5ad_dir <- dirname(output_h5ad)
dir.create(h5ad_dir, recursive = TRUE, showWarnings = FALSE)

ts("Writing .h5ad...")
adata$write_h5ad(output_h5ad)
ts("H5AD written to:", output_h5ad)
ts("All done.")
