message()

message('------------------------------------')
message("Starting: loading libraries...")
message('------------------------------------')

# ======================================================================
# == load libraries ==
# ======================================================================

suppressPackageStartupMessages({
  library(scater)
  library(Seurat)
  library(tidyverse)
  library(robustbase)
  library(patchwork)
  library(scDblFinder)
  library(presto)
  library(harmony)
})

# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# ======================================================================
# == parameters ==
# ======================================================================

# read the paramters from the config file
params    = snakemake@params$seurat_params
in_id_org = snakemake@params$organism_id
df_cc     = snakemake@params$df_cc

message('------------------------------------')
message("Organism Id:", in_id_org)
message('------------------------------------')

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# read in the rds file
rds               <- snakemake@input$rds
LUT_QC            <- snakemake@input$LUT_QC
# samples_meta_file <- snakemake@input$samples_meta

# define the output
out_merged_object     <- snakemake@output$merged_rds
out_integrated_object <- snakemake@output$integrated_rds
out_meta              <- snakemake@output$meta
out_markers           <- snakemake@output$markers

message('------------------------------------')
message("input files:"); 
message('------------------------------------')
cat(rds, sep = '\n')
message('------------------------------------')
message("output files:"); 
message('------------------------------------')
message("output merged:     ", out_merged_object)
message("output integrated: ", out_integrated_object)
message("output meta:       ", out_meta)
message("output markers:    ", out_markers)

# ======================================================================
# == custom functions ==
# ======================================================================

# run:
# - FindNeighbors
# - RunUMAP
# - FindClusters
# taking as input:
# - the reduction layer ('pca', 'harmony')
# - the parameters from snakemake

seurat_embedding = function(x, reduction, params) {
  x %>% 
    FindNeighbors(
      reduction = reduction, 
      # ======================================= #
      dims = 1:params$FindNeighbors$dims, 
      # ======================================= #
      assay = NULL,
      features = NULL, # Features to use as input for building the (S)NN; used only when dims is NULL
      # ======================================= #
      k.param = params$FindNeighbors$k.param, # Defines k for the k-nearest neighbor algorithm (SdP: change to make it compliant with RunUMAP)
      # ======================================= #
      return.neighbor = FALSE,
      # compute.SNN = !return.neighbor,
      prune.SNN = 1/15, # Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph
      nn.method = "annoy", # Method for nearest neighbor finding. Options include: rann, annoy
      n.trees = 50, # More trees gives higher precision when using annoy approximate nearest neighbor search
      annoy.metric = "euclidean", # Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
      nn.eps = 0, # Error bound when performing nearest neighbor seach using RANN; default of 0.0 implies exact nearest neighbor search
      verbose = TRUE,
      do.plot = FALSE,
      graph.name = NULL, # Default is assay.name_(s)nn
      l2.norm = FALSE, # Take L2Norm of the data
      cache.index = FALSE
    ) %>% 
    RunUMAP(
      # ======================================= #
      dims = 1:params$RunUMAP$dims, 
      # ======================================= #
      reduction = reduction,
      # ======================================= #
      features = NULL, # dims must be NULL to run on features
      graph = NULL, # Name of graph on which to run UMAP
      assay = DefaultAssay(object = object),
      nn.name = NULL, # Name of knn output on which to run UMAP
      slot = "data",
      umap.method = "uwot",
      reduction.model = NULL, # DimReduc object that contains the umap model
      return.model = FALSE,
      # ======================================= #
      n.neighbors = params$RunUMAP$n.neighbors, # SdP: same as in FindNeighbors?
      # ======================================= #
      n.components = 2L,
      metric = "cosine",
      n.epochs = NULL,
      # ======================================= #
      learning.rate = params$RunUMAP$learning.rate,
      # ======================================= #
      min.dist = params$RunUMAP$min.dist,
      # ======================================= #
      spread = params$RunUMAP$spread,
      # ======================================= #
      set.op.mix.ratio = 1,
      local.connectivity = 1L,
      repulsion.strength = 1,
      negative.sample.rate = 5L,
      a = NULL,
      b = NULL,
      uwot.sgd = FALSE,
      seed.use = 42L,
      metric.kwds = NULL,
      angular.rp.forest = FALSE,
      densmap = FALSE,
      dens.lambda = 2,
      dens.frac = 0.3,
      dens.var.shift = 0.1,
      verbose = TRUE,
      reduction.name = "umap",
      reduction.key = NULL
    ) %>% 
    FindClusters(
      graph.name = NULL,
      cluster.name = NULL,
      modularity.fxn = 1,
      initial.membership = NULL,
      node.sizes = NULL,
      # ======================================= #
      resolution = params$FindClusters$resolution, 
      # ======================================= #
      # argument is now depracated
      # method = "matrix",
      algorithm = 1,
      n.start = 10,
      n.iter = 10,
      random.seed = 0,
      group.singletons = TRUE,
      temp.file.location = NULL,
      edge.file.name = NULL,
      verbose = TRUE
    )
}

# ======================================================================
# == load and rename the barcodes ==
# ======================================================================

# message('------------------------------------')
# message("## Loading samples metadata...")
# message('------------------------------------')
# message()
# 
# samples_meta = read.csv(samples_meta_file)
# id_col = colnames(samples_meta)=='sample_name'
# annot_col = !id_col
# # rename the id colname
# colnames(samples_meta)[id_col] = 'orig.ident'
# # show in the log the annotation columns
# message('columnn: (values)')
# for( k in colnames(samples_meta)[annot_col] ) {
#   message(k, ' : ', paste(samples_meta[k], collapse = ','))
# }

message()
message('------------------------------------')
message("## Loading samples...")
message('------------------------------------')

# load objects
all_objects = lapply(rds, readRDS)

# ======================================================================
# == excluding samples flagged in the LUT_QC ==
# ======================================================================

# drop samples the user has flagged include_sample = FALSE in the (hand
# edited) LUT_QC table, before they ever reach merge()/integration. runApplyQC
# still ran the full per-sample QC pipeline for them, so their outputs remain
# available on disk for review even though they are excluded here.
df_lut <- read_csv(LUT_QC)
included_samples <- df_lut$sample_name[df_lut$include_sample]

sample_ids <- sapply(all_objects, Project)
excluded_samples <- setdiff(sample_ids, included_samples)
if (length(excluded_samples) > 0) {
  message("Excluding samples from merge/integration (include_sample = FALSE in LUT_QC): ",
          paste(excluded_samples, collapse = ", "))
}
all_objects <- all_objects[sample_ids %in% included_samples]

if (length(all_objects) == 0) {
  stop("All samples excluded via include_sample in LUT_QC -- nothing left to integrate.")
}

# rename cells (to avoid duplicated barcodes)
# fix the metadata to avoid confusion with resolutions names
lapply(all_objects, function(obj) {
  obj = RenameCells(obj, add.cell.id = Project(obj))
  
  # remove the clustering information to avoid confusions after integration
  meta_in <- obj@meta.data
  meta_fix <- meta_in %>%
    select(-contains("RNA_snn_res")) %>%
    select(-"seurat_clusters")
  
  # swap the metadata
  obj@meta.data <- meta_fix
  
  return(obj)
}) -> all_objects

# ======================================================================
# == filtering genes by sample  ==
# ======================================================================

message('------------------------------------')
message("## Filtering genes in samples...")
message('------------------------------------')
message('paramters:')
str(params['FilterGenes'])
message('------------------------------------')
min_fraction_of_cells = 
  params$FilterGenes$min_fraction_of_cells
all_objects = 
  lapply(all_objects, function(obj) {
    cells_per_gene = rowSums(GetAssayData(obj, layer = 'counts') > 0)
    cells_fraction_per_gene = cells_per_gene / ncol(obj)
    obj[cells_fraction_per_gene > min_fraction_of_cells,]
  })

# ======================================================================
# == displaying dataset features  ==
# ======================================================================

# count the number of objects to be integrated
n_objects = length(all_objects)
n_cells_bySample = sapply(all_objects,ncol)
n_genes_bySample = sapply(all_objects,nrow)
n_cells = sum(n_cells_bySample)
n_genes = length(unique(unlist(lapply(all_objects,rownames))))

message('------------------------------------')
message("Number of samples:\t", n_objects)
message("Number of cells:\t", n_cells, 
        ' (', paste(n_cells_bySample, collapse = ', '), ')')
message("Number of genes:\t", n_genes, 
        ' (', paste(n_genes_bySample, collapse = ', '), ')')
message('------------------------------------')

# ======================================================================
# == merge  ==
# ======================================================================

message('------------------------------------')
message("## Merging...")
message('------------------------------------')

# Ensures the limit is either your calculated size or 1GB, whichever is larger
dynamic_size <- n_objects * n_cells * 0.1 * 1024^2
options(future.globals.maxSize = max(dynamic_size, 1024^3))

# merge
merged = merge(all_objects[[1]], all_objects[-1])

# # -------- add samples meta to the merged object ------ #
# 
# # move rownames into a column
# single_cell_meta = merged@meta.data %>% 
#   tibble::rownames_to_column('barcode')
# # merge the two metadata
# single_cell_meta = 
#   merge(
#     single_cell_meta, 
#     samples_meta, 
#     by = 'orig.ident', 
#     all = T, sort = F) %>% 
#   tibble::column_to_rownames('barcode')
# # add the metadata annotations (be sure that the order is preserved)
# merged@meta.data = single_cell_meta[colnames(merged),]

# ======================================================================
# == delete all_objects and free up the memory ==
# ======================================================================

message('------------------------------------')
message("## Completed. Removing samples...")
message('------------------------------------')

rm(all_objects); gc()

# ======================================================================
# == preprocess merged data ==
# ======================================================================

message('------------------------------------')
message("## Normalizing merged object...")
message('------------------------------------')
message('paramters:')
str(params['NormalizeData'])
message('------------------------------------')
merged <- merged %>% 
  NormalizeData(
    assay = NULL,
    normalization.method = params$NormalizeData$normalization.method,
    scale.factor = params$NormalizeData$scale.factor,
    margin = 1,
    verbose = TRUE
  ) 

# ======================================================================
# == calculate CellCycle score ==
# ======================================================================

message('------------------------------------')
message("## Adding CellCycle scores...")
message('------------------------------------')

# choose the right set of genes
if(in_id_org == "9606"){
  cell.cycle.genes = read.table(
    df_cc[["9606"]], 
    sep = '\t', header = T)
} else if(in_id_org == "10090"){
  cell.cycle.genes = read.table(
    df_cc[["10090"]], 
    sep = '\t', header = T)
} else {
  stop('Organism specific cell-cycle genes not found.')
}

# add the mefadata accordin to the organism
s.genes = cell.cycle.genes %>% 
  dplyr::filter(class == 's.genes') %>% 
  pull(gene)
g2m.genes = cell.cycle.genes %>% 
  dplyr::filter(class == 'g2m.genes') %>% 
  pull(gene)
merged <- CellCycleScoring(merged, 
                           s.features = s.genes, 
                           g2m.features = g2m.genes)

# ======================================================================
# == Scaling and PCA ==
# ======================================================================

message('------------------------------------')
message("## Scaling and calculating PCs ...")
message('------------------------------------')
message('paramters:')
str(params[c('FindVariableFeatures','ScaleData','RunPCA')])
message('------------------------------------')
merged <- merged %>% 
  FindVariableFeatures(
    assay = NULL,
    selection.method = params$FindVariableFeatures$selection.method,
    loess.span = 0.3,
    clip.max = "auto",
    # mean.function = FastExpMean,
    # dispersion.function = FastLogVMR,
    num.bin = 20,
    binning.method = "equal_width",
    nfeatures = params$FindVariableFeatures$nfeatures,
    mean.cutoff = c(0.1, 8),
    dispersion.cutoff = c(1, Inf),
    verbose = TRUE
  ) %>%
  ScaleData(
    features = NULL, # Default is variable features.
    assay = NULL,
    vars.to.regress = params$ScaleData$vars.to.regress,
    split.by = NULL,
    model.use = "linear",
    use.umi = FALSE,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10,
    block.size = 1000,
    min.cells.to.block = 3000,
    verbose = FALSE
  ) %>%
  RunPCA(
    assay = NULL,
    features = NULL, # If features=NULL, PCA will be run using the variable features for the Assay
    npcs = params$RunPCA$npcs,
    rev.pca = FALSE,
    weight.by.var = TRUE,
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = 42
  ) 

# ======================================================================
# == Embedding ==
# ======================================================================

message('------------------------------------')
message("## Embedding...")
message('------------------------------------')
message('paramters:')
str(params[c('FindNeighbors','RunUMAP','FindClusters')])
# embed based on pca reduction
merged <- merged %>% 
  seurat_embedding(
    reduction = 'pca', 
    params = params
  )

# ======================================================================
# == integrate via Harmony ==
# ======================================================================

message('------------------------------------')
message("## Harmony integration...")
message('------------------------------------')
message('paramters:')
str(params['HarmonyIntegration'])
message('------------------------------------')
if( n_objects > 1 ) # need integration
{
  integrated = IntegrateLayers(
    merged, 
    method = HarmonyIntegration,
    orig.reduction = "pca",
    assay = NULL,
    features = NULL, # Ignored
    scale.layer = "scale.data", # Ignored
    new.reduction = "harmony",
    layers = NULL, # Ignored
    npcs = NULL,
    key = "harmony_",
    theta = params$HarmonyIntegration$theta, # Diversity clustering penalty parameter
    lambda = params$HarmonyIntegration$lambda, # Ridge regression penalty parameter
    sigma = params$HarmonyIntegration$sigma, # Width of soft kmeans clusters
    nclust = NULL, # Number of clusters in model
    tau = 0, # Protection against overclustering small datasets with large ones
    block.size = 0.05, # What proportion of cells to update during clustering
    max.iter.harmony = params$HarmonyIntegration$max.iter.harmony, # Maximum number of rounds to run Harmony
    max.iter.cluster = 20L, # Maximum number of rounds to run clustering at each round of Harmony
    epsilon.cluster = 1e-05, # Convergence tolerance for clustering round of Harmony
    epsilon.harmony = 1e-04, # Convergence tolerance for Harmony
    verbose = TRUE
  )
} else { # no-need for integration
  integrated = merged
  integrated@reductions$harmony = integrated@reductions$pca
}

# ======================================================================
# == save and erase merged object ==
# ======================================================================

message('------------------------------------')
message("## Integration completed.")
message('------------------------------------')

# join and remove stale layers of scale.data
merged = JoinLayers(merged)
stale_layers <- grep("^scale\\.data\\.", Layers(merged), value = TRUE)
for (lyr in stale_layers) {
  LayerData(merged, assay = "RNA", layer = lyr) <- NULL
}

message('------------------------------------')
message("## Saving and removing merged...")
message('------------------------------------')

saveRDS(merged, out_merged_object)
rm(merged); gc()

# ======================================================================
# == embed and cluster ==
# ======================================================================

message('------------------------------------')
message("## Embedding harmonized object...")
message('------------------------------------')
message('paramters:')
str(params[c('FindNeighbors','RunUMAP','FindClusters')])
message('------------------------------------')
# embed based on hermony reduction
integrated <- integrated %>% 
  seurat_embedding(
    reduction = "harmony", 
    params = params
  )

# ======================================================================
# == save integrated object ==
# ======================================================================

message('------------------------------------')
message("## Saving harmonized object...")
message('------------------------------------')

integrated = JoinLayers(integrated)
stale_layers <- grep("^scale\\.data\\.", Layers(integrated), value = TRUE)
for (lyr in stale_layers) {
  LayerData(integrated, assay = "RNA", layer = lyr) <- NULL
}
integrated@misc$seurat_params = params
saveRDS(integrated, out_integrated_object)
write_tsv(integrated@meta.data %>% 
            rownames_to_column("barcodes"), out_meta)

# ======================================================================
# == markers ==
# ======================================================================

message('------------------------------------')
message("## Calculating markers...")
message('------------------------------------')
message('paramters:')
str(params['FindAllMarkers'])
message('------------------------------------')
integrated.markers = 
  lapply(params$FindClusters$resolution, 
         function(res) {
           message('------- resolution: ', res)
           # set identity to the desired clustering resolution
           Idents(integrated) = integrated[[paste0(DefaultAssay(integrated), '_snn_res.',res)]][[1]]
           markers.res = FindAllMarkers(
             integrated, 
             assay = NULL,
             features = NULL,
             logfc.threshold = params$FindAllMarkers$logfc.threshold, # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
             test.use = params$FindAllMarkers$test.use,
             slot = "data",
             min.pct = params$FindAllMarkers$min.pct, # SdP: default is 0.01
             min.diff.pct = -Inf, # only test genes that show a minimum difference in the fraction of detection between the two groups.
             node = NULL, # A node to find markers for and all its children; requires BuildClusterTree to have been run previously; replaces FindAllMarkersNode
             verbose = TRUE,
             only.pos = params$FindAllMarkers$only.pos, # SdP: default is FALSE
             max.cells.per.ident = Inf, # Down sample each identity class to a max number. Default is no downsampling
             random.seed = 1, # Random seed for downsampling
             latent.vars = NULL, # Variables to test, used only when test.use is one of 'LR', 'negbinom', 'poisson', or 'MAST'
             min.cells.feature = 3, # Minimum number of cells expressing the feature in at least one of the two groups, currently only used for poisson and negative binomial tests
             min.cells.group = 3, # Minimum number of cells in one of the groups
             mean.fxn = NULL, # Function to use for fold change or average difference calculation. The default depends on the the value of fc.slot
             fc.name = NULL, # Name of the fold change, average difference, or custom function column in the output data.frame. If NULL, the fold change column will be named according to the logarithm base (eg, "avg_log2FC"), or if using the scale.data slot "avg_diff".
             base = 2, # The base with respect to which logarithms are computed.
             return.thresh = 0.01, # Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
             densify = FALSE # Convert the sparse matrix to a dense form before running the DE test. This can provide speedups but might require higher memory; default is FALSE   
           )
           markers.res %>% 
             mutate(diff.pct = pct.1 - pct.2) %>%
             mutate(resolution = res)
         })
# row binding to obtain a unique data.frame
integrated.markers = do.call('rbind', integrated.markers)

# ======================================================================
# == write markers ==
# ======================================================================

message('------------------------------------')
message("## Saving markers...")
message('------------------------------------')

write_tsv(integrated.markers, out_markers)

message('------------------------------------')
message("## Done.")
message('------------------------------------')
