# ======================================================================
# == load libraries ==
# ======================================================================

message('------------------------------------')
message("Starting: loading libraries...")
message('------------------------------------')

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

# ======================================================================
# == custom functions ==
# ======================================================================

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# ======================================================================
# == parameters ==
# ======================================================================

# inputs
integrated_path  = snakemake@input$rds
annotation_paths = snakemake@input$annotation_metas

# outputs
integrated_annot_path   = snakemake@output$rds

# ======================================================================
# == load integrated ==
# ======================================================================

message('------------------------------------')
message("Loading the integrated dataset...")
message('------------------------------------')

integrated = readRDS(integrated_path)

# ======================================================================
# == add the annotation(s) ==
# ======================================================================

message('------------------------------------')
message("Adding annotations...")
message('------------------------------------')

integrated = reduce(annotation_paths, function(obj, annotation_path) {
  message("  merging: ", annotation_path)
  annotation_meta <- read.csv(annotation_path) %>%
    column_to_rownames("X")
  AddMetaData(obj, annotation_meta)
}, .init = integrated)

# ======================================================================
# == save annotated integrated ==
# ======================================================================

message('------------------------------------')
message("Saving the annotated integrated dataset...")
message('------------------------------------')

saveRDS(integrated, integrated_annot_path)
