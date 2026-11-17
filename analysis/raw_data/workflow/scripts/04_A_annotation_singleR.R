# ======================================================================
# == load libraries ==
# ======================================================================

message('------------------------------------')
message("Starting: loading libraries...")
message('------------------------------------')

suppressPackageStartupMessages({
  library(scRNAseq)
  library(SingleR)
  library(Seurat)
  library(dplyr)
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
integrated_path   = snakemake@input$rds

# params
reference_dataset = snakemake@params$reference_dataset

# outputs
annotation_meta   = snakemake@output$annotation_meta
singleR_output    = c(
  'main' = snakemake@output$singleR_output_main,
  'fine' = snakemake@output$singleR_output_fine
)
singleR_heatmap   = c(
  'main' = snakemake@output$singleR_heatmap_main,
  'fine' = snakemake@output$singleR_heatmap_fine
)

# ======================================================================
# == load integrated ==
# ======================================================================

message('------------------------------------')
message("Loading the integrated dataset...")
message('------------------------------------')

integrated = readRDS(integrated_path)

# ======================================================================
# == singleR reference ==
# ======================================================================

message('------------------------------------')
message("Loading SingleR reference...")
message('------------------------------------')

ref_fun <- get(reference_dataset, envir = asNamespace("celldex"))
ref <- ref_fun()

# lognormalize if necessary
if(! "logcounts" %in% names(assays(ref))) {
  ref = scuttle::logNormCounts(ref)
}

# ======================================================================
# == make prediction ==
# ======================================================================

prediction_types = c('fine', 'main')

pred.integrated.list = 
  
  lapply(prediction_types, function(label) {

    message('------------------------------------')
    message(paste("SingleR", label, "prediction..."))
    message('------------------------------------')
    
    pred.integrated <- SingleR(
      test = GetAssayData(
        integrated, 
        assay = DefaultAssay(integrated), 
        layer = 'data'), 
      ref = ref, 
      assay.type.ref='logcounts',
      labels = ref@colData[[paste("label", label, sep = '.')]]
      )
    
    saveRDS(pred.integrated, file = singleR_output[label])

    plotScoreHeatmap(pred.integrated,
                     filename = singleR_heatmap[label],
                     height=6, width = 10)
    
    return(pred.integrated)
    
  })

message('------------------------------------')
message("Writing prediction to file...")
message('------------------------------------')

# create a data.frame with the predictions
class_list      = sapply(pred.integrated.list, '[', 'labels')
class_dataframe = merge(class_list[[1]], class_list[[2]], by = 0) %>% 
  data.frame %>% 
  tibble::column_to_rownames('Row.names')
colnames(class_dataframe) = paste('CELLTYPE_SingleR', prediction_types, sep = '_')

# create a data.frame with the quality of the predictions
qual_list      = sapply(pred.integrated.list, '[', 'delta.next')
qual_dataframe = merge(qual_list[[1]], qual_list[[2]], by = 0) %>% 
  data.frame %>% 
  tibble::column_to_rownames('Row.names')
colnames(qual_dataframe) = paste('CELLTYPEQ_SingleR', prediction_types, sep = '_')

# combine the 2
pred_dataframe = cbind(class_dataframe, qual_dataframe)

# write to file
write.csv(pred_dataframe, annotation_meta)
