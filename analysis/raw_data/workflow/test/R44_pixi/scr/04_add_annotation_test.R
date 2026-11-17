# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# parameters --------------------------------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# inputs
# integrated_path   = snakemake@input$rds
integrated_path <- "data/integrated_obj_default.rds"
# annotation_paths = snakemake@input$annotation_metas
annotation_paths <- c(
  singleR   = "out/table/integrated_singleR_meta_default.csv",
  cytetypeR = "out/table/integrated_cytetypeR_meta_default.csv"
)

# outputs
# integrated_annot_path   = snakemake@output$rds
integrated_annot_path <- "out/object/integrated_obj_default_annotated.rds"

# load integrated ---------------------------------------------------------
integrated <- readRDS(integrated_path)

# add annotation(s) ---------------------------------------------------------
message('------------------------------------')
message("Adding annotations...")
message('------------------------------------')

# fail fast if any table is missing, before merging any of them
walk(annotation_paths, function(annotation_path) {
  if (!file.exists(annotation_path)) {
    stop("Annotation table not found: ", annotation_path,
         " -- run the corresponding annotation test script first.")
  }
})

integrated <- purrr::reduce(annotation_paths, .init = integrated, function(obj, annotation_path) {
  message("  merging: ", annotation_path)
  annotation_meta <- read.csv(annotation_path) %>%
    column_to_rownames("X")
  AddMetaData(obj, annotation_meta)
})

integrated@meta.data %>% head()

# ======================================================================
# == check the effect of integrating both annotation methods ==
# ======================================================================

message('------------------------------------')
message("Inspecting merged annotations...")
message('------------------------------------')

singleR_cols   <- grep("^CELLTYPEQ_SingleR",  colnames(integrated@meta.data), value = TRUE)
cytetypeR_cols <- grep("^cytetype_",         colnames(integrated@meta.data), value = TRUE)

message("SingleR columns added: ",   paste(singleR_cols,   collapse = ", "))
message("CyteTypeR columns added: ", paste(cytetypeR_cols, collapse = ", "))

# every cell should have been assigned by both methods, since both run on
# the full integrated object -- flag it loudly if that's not the case
message("Cells missing SingleR annotation: ",
        sum(is.na(integrated@meta.data[[singleR_cols[1]]])))
message("Cells missing CyteTypeR annotation: ",
        sum(is.na(integrated@meta.data[[cytetypeR_cols[1]]])))

# cross-tabulate the two "main" cell-type calls to see how much they agree
cytetypeR_ann_col <- grep("^cytetype_annotation", cytetypeR_cols, value = TRUE)[1]
print(table(
  SingleR   = integrated@meta.data[["CELLTYPE_SingleR_main"]],
  CyteTypeR = integrated@meta.data[[cytetypeR_ann_col]]
))

# ======================================================================
# == save annotated integrated ==
# ======================================================================

message('------------------------------------')
message("Saving the annotated integrated dataset...")
message('------------------------------------')

saveRDS(integrated, integrated_annot_path)
