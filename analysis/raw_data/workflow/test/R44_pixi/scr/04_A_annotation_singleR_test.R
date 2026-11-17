# libraries ---------------------------------------------------------------
library(scRNAseq)
library(SingleR)
library(Seurat)
library(tidyverse)
library(celldex)

# - r-seurat=5.3.1
# - r-tidyverse=2.0.0
# - r-robustbase=0.99_6
# - bioconductor-scuttle=1.16.0
# - bioconda::bioconductor-singler=2.8.0
# - bioconda::bioconductor-scrnaseq=2.20.0
# - bioconda::bioconductor-celldex=1.16.0
# - conda-forge::r-viridis=0.6.5
# - conda-forge::r-pheatmap=1.0.13

# parameters --------------------------------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# inputs
# integrated_path   = snakemake@input$rds
integrated_path <- "data/integrated_obj_default.rds"

# params
# reference_dataset = snakemake@params$reference_dataset
reference_dataset <- "BlueprintEncodeData"

# outputs
# annotation_meta <- snakemake@output$annotation_meta
annotation_meta <- "out/table/integrated_singleR_meta_default.csv"
# singleR_output <- c(
#   'main' = snakemake@output$singleR_output_main,
#   'fine' = snakemake@output$singleR_output_fine)
# singleR_heatmap <- c(
#   'main' = snakemake@output$singleR_heatmap_main,
#   'fine' = snakemake@output$singleR_heatmap_fine)

singleR_output <- c(
  'main' = "out/object/integrated_singleR_main_prediction_default.rds",
  'fine' = "out/object/integrated_singleR_fine_prediction_default.rds")

singleR_heatmap <- c(
  'main' = "out/plot/integrated_singleR_main_plotScoreHeatmap_default.png",
  'fine' = "out/plot/integrated_singleR_fine_plotScoreHeatmap_default.png")

# load integrated ---------------------------------------------------------
integrated <- readRDS(integrated_path)
integrated@meta.data

DimPlot(integrated)

# singleR reference -------------------------------------------------------
ref_fun <- get(reference_dataset, envir = asNamespace("celldex"))
ref <- ref_fun()

# lognormalize if necessary
if(! "logcounts" %in% names(assays(ref))) {
  ref = scuttle::logNormCounts(ref)
}

# make prediction ---------------------------------------------------------
prediction_types = c('fine', 'main')

# label <- "main"
pred.integrated.list <- 
  
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
    
    # save annotation barcodes
    saveRDS(pred.integrated, file = singleR_output[label])
    
    # save the plot
    plotScoreHeatmap(pred.integrated,
                     filename = singleR_heatmap[label],
                     height=6, width = 10)
    
    return(pred.integrated)
    
  })

message('------------------------------------')
message("Writing prediction to file...")
message('------------------------------------')

# create a data.frame with the predictions
class_list <- sapply(pred.integrated.list, '[', 'labels')
class_dataframe <- merge(class_list[[1]], class_list[[2]], by = 0) %>% 
  data.frame %>% 
  tibble::column_to_rownames('Row.names')
colnames(class_dataframe) <- paste('CELLTYPE_SingleR', prediction_types, sep = '_')

# create a data.frame with the quality of the predictions
qual_list <- sapply(pred.integrated.list, '[', 'delta.next')
qual_dataframe <- merge(qual_list[[1]], qual_list[[2]], by = 0) %>% 
  data.frame %>% 
  tibble::column_to_rownames('Row.names')
colnames(qual_dataframe) = paste('CELLTYPEQ_SingleR', prediction_types, sep = '_')

# combine the 2
pred_dataframe <- cbind(class_dataframe, qual_dataframe)

# write to file
write.csv(pred_dataframe, annotation_meta)
