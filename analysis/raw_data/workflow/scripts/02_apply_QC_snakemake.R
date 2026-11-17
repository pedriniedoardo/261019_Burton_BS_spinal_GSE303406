# ======================================================================
# == load libraries ==
# ======================================================================

library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(scDblFinder)

# ======================================================================
# == custom functions ==
# ======================================================================

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# read in the rds file
rds <- snakemake@input$rds

# sample id
id_sample <- snakemake@wildcards$sample_name

# read in the LUT_QC
LUT_QC <- snakemake@input$LUT_QC

message("input rds: ", rds)
message("sample id: ", id_sample)
message("table LUT QC: ", LUT_QC)

# define the output
# out_id_object <-  "results/Seurat/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
# out_id_meta <- "results/Seurat/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.rds"
out_id_object <- snakemake@output$rds
out_id_meta <- snakemake@output$meta

message("output rds: ", out_id_object)
message("output meta: ", out_id_meta)


# ======================================================================
# == load input ==
# ======================================================================

df_LUT <- read_csv(LUT_QC) %>%
  filter(sample_name == id_sample)

nCount_RNA_high <- df_LUT$nCount_RNA_high
nCount_RNA_low <- df_LUT$nCount_RNA_low

nFeature_RNA_high <- df_LUT$nFeature_RNA_high
nFeature_RNA_low <- df_LUT$nFeature_RNA_low

percent.ribo_high <- df_LUT$percent.ribo_high
percent.ribo_low <- df_LUT$percent.ribo_low

percent.mt_high <- df_LUT$percent.mt_high
percent.mt_low <- df_LUT$percent.mt_low

percent.globin_high <- df_LUT$percent.globin_high
percent.globin_low <- df_LUT$percent.globin_low

# pull the parameter for filtering the doublets
remove_doublets <- df_LUT$remove_doublets

# build the label
label <- paste("nCount-",nCount_RNA_low,"_",nCount_RNA_high,"|",
               "nFeat-",nFeature_RNA_low,"_",nFeature_RNA_high,"|",
               "pRibo-",percent.ribo_low,"_",percent.ribo_high,"|",
               "pMito-",percent.mt_low,"_",percent.mt_high,"|",
               "pHemo-",percent.globin_low,"_",percent.globin_high,"|",
               "wDbl-",remove_doublets,"|",
               "V5",sep = "")

# ======================================================================
# == run the standard processing ==
# ======================================================================

scobj <- readRDS(rds)

# add the filtering label
scobj$label <- label

# add the filtering variable based on the fixed threshold
scobj$discard_threshold <- scobj@meta.data %>%
  mutate(test = nCount_RNA < nCount_RNA_low | nCount_RNA > nCount_RNA_high |
                nFeature_RNA < nFeature_RNA_low | nFeature_RNA > nFeature_RNA_high |
                percent.ribo < percent.ribo_low | percent.ribo > percent.ribo_high |
                percent.mt < percent.mt_low | percent.mt > percent.mt_high |
                percent.globin < percent.globin_low | percent.globin > percent.globin_high) %>%
  pull(test)

# preprocess the dataset before the doublet identification as recommended in:
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
# 1.5.11
# according to the documentation the doublet identification should be run before the fine filtering
# remove low coverage cells and proprocess to generata clusters needed for the doublet identification
scobj <- subset(scobj,subset = nCount_RNA > 500) %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  # do I run the UMAP ? I do not need it for the doublet identification, but can be useful in case someone wants to explore an individual sample
  RunUMAP(dims = 1:30)

# run scDblFinder after filtering the low coverage cells
sce_scobj <- scDblFinder(GetAssayData(scobj, layer="counts"), clusters=Idents(scobj))

# port the resulting scores back to the Seurat object:
scobj$scDblFinder.score <- sce_scobj$scDblFinder.score
scobj$scDblFinder.class <- sce_scobj$scDblFinder.class

# cross-tabulate doublet calls against the threshold-based discard flag
table(scobj$scDblFinder.class,scobj$discard_threshold)

# perform the filtering based on the fixed threshold defined, and on doublets if requested
if(remove_doublets){
  # perform the filtering based on the fixed threshold defined and for the doublets
  scobj_filter <- subset(scobj, subset = discard_threshold == 0 & scDblFinder.class == "singlet")
}else{
  # perform the filtering based on the fixed threshold defined only
  scobj_filter <- subset(scobj, subset = discard_threshold == 0)
}

# preprocess data after filtering 
scobj_filter <- scobj_filter %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30) 

# ======================================================================
# == save output ==
# ======================================================================

saveRDS(scobj_filter,out_id_object)
write_tsv(scobj_filter@meta.data %>% rownames_to_column("barcodes"),out_id_meta)