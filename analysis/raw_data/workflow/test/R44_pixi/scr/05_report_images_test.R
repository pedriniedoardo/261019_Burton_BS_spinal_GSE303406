# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)

# parameters --------------------------------------------------------------
options(Seurat.object.assay.version = "v5")

# focused test of the annotation-column discovery change in
# 05_report_images.R (celltype_meta grep): must pick up SingleR's
# CELLTYPE_ columns and CyteTypeR's annotation/cellOntologyTerm/cellState
# columns, while excluding cellOntologyTermID and granularAnnotation --
# and must adapt automatically to whichever method(s) were run
# (config["annotation_methods"] switch happens upstream, in addAnnotations).

integrated_paths <- c(
  both    = "data/integrated_obj_cellbender_annotated.rds",
  onlySR  = "data/integrated_obj_cellbender_annotated_onlySR.rds",
  onlyCT  = "data/integrated_obj_cellbender_annotated_onlyCT.rds",
  none  = "data/integrated_obj_cellbender_annotated_none.rds"
)

out_location <- "out/"
figure_extension <- "png"

# the exact pattern used in 05_report_images.R --------------------------
celltype_pattern <- 'CELLTYPE_|cytetype_annotation_|cytetype_cellOntologyTerm(?!ID)|cytetype_cellState_'

# scenario <- "both"
walk(names(integrated_paths), function(scenario) {

  message('------------------------------------')
  message("Scenario: ", scenario, " (", integrated_paths[scenario], ")")
  message('------------------------------------')

  integrated <- readRDS(integrated_paths[scenario])

  all_ct_related <- grep('CELLTYPE|cytetype', colnames(integrated@meta.data),
                          value = T, ignore.case = T)
  celltype_meta <- grep(celltype_pattern, colnames(integrated@meta.data),
                         value = T, perl = T)
  celltypeq_meta <- grep('CELLTYPEQ_', colnames(integrated@meta.data),
                          value = T)

  message("All CELLTYPE/cytetype columns present: ",
          paste(all_ct_related, collapse = ", "))
  message("Selected for categorical DimPlot (celltype_meta): ",
          paste(celltype_meta, collapse = ", "))
  message("Selected for numeric FeaturePlot (celltypeq_meta): ",
          paste(celltypeq_meta, collapse = ", "))

  # sanity checks -----------------------------------------------------
  excluded <- setdiff(all_ct_related, c(celltype_meta, celltypeq_meta))
  message("Columns intentionally excluded: ", paste(excluded, collapse = ", "))
  stopifnot(!any(grepl("cellOntologyTermID", celltype_meta)))
  stopifnot(!any(grepl("granularAnnotation", celltype_meta)))

  # render one DimPlot per selected column, to confirm they actually plot
  for (ct_meta in celltype_meta) {
    p <- DimPlot(integrated, group.by = ct_meta, label = TRUE) +
      ggtitle(paste(scenario, ct_meta, sep = " | "))
    ggsave(
      paste0(out_location, "plot/dim_test_", scenario, "_", ct_meta, ".", figure_extension),
      p, height = 7, width = 9
    )
  }

  for (ctq_meta in celltypeq_meta) {
    p <- FeaturePlot(integrated, features = ctq_meta) +
      ggtitle(paste(scenario, ctq_meta, sep = " | "))
    ggsave(
      paste0(out_location, "plot/dim_test_", scenario, "_", ctq_meta, ".", figure_extension),
      p, height = 6, width = 8
    )
  }

  message("Plots written for scenario '", scenario, "'.")
  message()
})

message('------------------------------------')
message("Done. Check ", out_location, "plot/ for dim_test_* files per scenario.")
message('------------------------------------')
