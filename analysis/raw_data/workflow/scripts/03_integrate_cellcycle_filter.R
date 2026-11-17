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

# ======================================================================
# == custom functions ==
# ======================================================================

# library(homologene)
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# s.genes <- homologene(s.genes, inTax = 9606, outTax = 10090)$`10090`
# g2m.genes <- homologene(g2m.genes, inTax = 9606, outTax = 10090)$`10090`
# 
# bind_rows(data.frame(gene = c(cc.genes$s.genes,cc.genes$g2m.genes),
#                      class = c(rep("s.genes",length(cc.genes$s.genes)),
#                                rep("g2m.genes",length(cc.genes$g2m.genes)))) %>%
#             mutate(org = "9606"),
#           data.frame(gene = c(s.genes,g2m.genes),
#                      class = c(rep("s.genes",length(s.genes)),
#                                rep("g2m.genes",length(g2m.genes)))) %>%
#             mutate(org = "10090")) %>%
#   write_tsv("../../../metadata/data.frame_cc.tsv")

# data_cc = "../../../metadata/data.frame_cc.tsv"
# org_id = "9606"
score_cell_cycle <- function(seurat_obj, org_id = "9606", data_cc = "../../../metadata/data.frame_cc.tsv") {
  
  if (org_id == "9606") {
    message("processing human genes")
    # load the cc genes for the human dataset
    df <- read_tsv(data_cc) %>%
      filter(org == org_id)
    s.genes <- df %>%
      filter(class == "s.genes") %>%
      pull(gene) %>% 
      unique()
    g2m.genes <- df %>%
      filter(class == "g2m.genes") %>%
      pull(gene) %>% 
      unique()
  } else if (org_id == "10090") {
    message("processing mouse genes")
    # read in the cc genes for the mouse dataset
    df <- read_tsv(data_cc) %>%
      filter(org == org_id)
    s.genes <- df %>%
      filter(class == "s.genes") %>%
      pull(gene) %>% 
      unique()
    g2m.genes <- df %>%
      filter(class == "g2m.genes") %>%
      pull(gene) %>% 
      unique()
  } else if (org_id != "9606") {
    stop("Species must be 'human' or 'mouse'")
  }
  
  # Run Seurat scoring
  seurat_obj <- CellCycleScoring(seurat_obj, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes)
  return(seurat_obj)
}

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# ======================================================================
# == parameters ==
# ======================================================================

n_dims = 30
resolutions = c(.3, .6, 1, 1.2)
min_pct = 0.1
cut_off_genes <- 20

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# read in the rds file
# rds <- c("data/sample_obj_cossu_default/object/Unsorted_PAO12790_default_obj_postQC.rds", "data/sample_obj_cossu_default/object/CD56_Plus_PAO12790_default_obj_postQC.rds")
rds <- snakemake@input$rds
# in_data_cc = "../../../metadata/data.frame_cc.tsv"
in_data_cc <- snakemake@input$df_cc

# define the output
out_merged_object     <- snakemake@output$merged_rds
out_integrated_object <- snakemake@output$integrated_rds
out_meta              <- snakemake@output$meta
out_markers           <- snakemake@output$markers

# define ogganisim
# in_id_org <- "9606"
in_id_org <- snakemake@params$id_org

message("input rds:"); cat(rds, sep = '\n')
message("input df_cc: ",in_data_cc)
message("output merged rds: ", out_merged_object)
message("output integrated rds: ", out_integrated_object)
message("output meta: ", out_meta)
message("output markers: ", out_markers)
message("id organism: ", in_id_org)

# ======================================================================
# == load and rename the barcodes ==
# ======================================================================

# load objects 
all_objects = lapply(rds, readRDS)

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
# == merge and pre-process ==
# ======================================================================

n_cells = sum(sapply(all_objects,ncol))

# Ensures the limit is either your calculated size or 1GB, whichever is larger
dynamic_size <- length(all_objects) * n_cells * 0.1 * 1024^2
options(future.globals.maxSize = max(dynamic_size, 1024^3))

# merge
# merged = merge(.all_objects[[1]], .all_objects[-1])

# test filtering low expressed genes --------------------------------------
# merge
merged_full <- merge(all_objects[[1]], all_objects[-1])
# apply the gene filtering after merging
# manual implementation
counts <- GetAssayData(JoinLayers(merged_full[["RNA"]]), assay = "RNA", layer = "counts")
genes_to_keep <- rowSums(counts > 0) >= cut_off_genes
merged <- merged_full[genes_to_keep, ]

# remove the object that are not needed anymore
rm(merged_full,counts,genes_to_keep)
gc()
# -------------------------------------------------------------------------

# standard pre-processing
# (-defaults- arguments exposed)
merged <- merged %>% 
  NormalizeData(
    assay = NULL,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    margin = 1,
    verbose = TRUE
  )

# add the cell cycle scoring
merged <- score_cell_cycle(merged,
                           org_id = in_id_org,
                           data_cc = in_data_cc)

# run the rest of the analysis
merged <- merged %>% 
  FindVariableFeatures(
    assay = NULL,
    selection.method = "vst",
    loess.span = 0.3,
    clip.max = "auto",
    # mean.function = FastExpMean,
    # dispersion.function = FastLogVMR,
    num.bin = 20,
    binning.method = "equal_width",
    nfeatures = 2000,
    mean.cutoff = c(0.1, 8),
    dispersion.cutoff = c(1, Inf),
    verbose = TRUE
  ) %>%
  ScaleData(
    features = NULL, # Default is variable features.
    assay = NULL,
    # vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"),
    vars.to.regress = NULL,
    split.by = NULL,
    model.use = "linear",
    use.umi = FALSE,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10,
    block.size = 1000,
    min.cells.to.block = 3000,
    verbose = TRUE
  ) %>%
  RunPCA(
    assay = NULL,
    features = NULL, # If features=NULL, PCA will be run using the variable features for the Assay
    npcs = 50,
    rev.pca = FALSE,
    weight.by.var = TRUE,
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = 42
  ) %>% 
  FindNeighbors(
    reduction = "pca", # SdP: cannot be changed
    dims = 1:n_dims, 
    assay = NULL,
    features = NULL, # Features to use as input for building the (S)NN; used only when dims is NULL
    k.param = 20, # Defines k for the k-nearest neighbor algorithm (SdP: change to make it compliant with RunUMAP)
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
    dims = 1:n_dims, 
    reduction = "pca", # SdP: cannot be changed
    features = NULL, # dims must be NULL to run on features
    graph = NULL, # Name of graph on which to run UMAP
    assay = DefaultAssay(object = object),
    nn.name = NULL, # Name of knn output on which to run UMAP
    slot = "data",
    umap.method = "uwot",
    reduction.model = NULL, # DimReduc object that contains the umap model
    return.model = FALSE,
    n.neighbors = 30L, # SdP: same as in FindNeighbors?
    n.components = 2L,
    metric = "cosine",
    n.epochs = NULL,
    learning.rate = 1,
    min.dist = 0.3,
    spread = 1,
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
    resolution = resolutions, # SdP: multiple resolutions?
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

# ======================================================================
# == integrate via Harmony ==
# ======================================================================

if( length(all_objects)>1 ) # need integration
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
    npcs = 50L,
    key = "harmony_",
    theta = NULL, # Diversity clustering penalty parameter
    lambda = NULL, # Ridge regression penalty parameter
    sigma = 0.1, # Width of soft kmeans clusters
    nclust = NULL, # Number of clusters in model
    tau = 0, # Protection against overclustering small datasets with large ones
    block.size = 0.05, # What proportion of cells to update during clustering
    max.iter.harmony = 10L, # Maximum number of rounds to run Harmony
    max.iter.cluster = 20L, # Maximum number of rounds to run clustering at each round of Harmony
    epsilon.cluster = 1e-05, # Convergence tolerance for clustering round of Harmony
    epsilon.harmony = 1e-04, # Convergence tolerance for Harmony
    verbose = TRUE
  )
} else { # no-need for integration
  integrated = merged
  integrated@reductions$harmony = integrated@reductions$pca
}

# remove the object that are not needed anymore
rm(all_objects)
gc()

# ======================================================================
# == embed and cluster ==
# ======================================================================

integrated <- integrated %>% 
  FindNeighbors(
    reduction = "harmony", # SdP: cannot be changed
    dims = 1:n_dims, 
    assay = NULL,
    features = NULL, # Features to use as input for building the (S)NN; used only when dims is NULL
    k.param = 20, # Defines k for the k-nearest neighbor algorithm (SdP: change to make it compliant with RunUMAP)
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
    dims = 1:n_dims, 
    reduction = "harmony", # SdP: cannot be changed
    features = NULL, # dims must be NULL to run on features
    graph = NULL, # Name of graph on which to run UMAP
    assay = DefaultAssay(object = object),
    nn.name = NULL, # Name of knn output on which to run UMAP
    slot = "data",
    umap.method = "uwot",
    reduction.model = NULL, # DimReduc object that contains the umap model
    return.model = FALSE,
    n.neighbors = 30L, # SdP: same as in FindNeighbors?
    n.components = 2L,
    metric = "cosine",
    n.epochs = NULL,
    learning.rate = 1,
    min.dist = 0.3,
    spread = 1,
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
    resolution = resolutions, # SdP: multiple resolutions?
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

# ======================================================================
# == join the layers of the merged object ==
# ======================================================================

merged     = JoinLayers(merged)
integrated = JoinLayers(integrated)

# remove stale layers -----------------------------------------------------

# format(object.size(integrated), units = "GB")
# format(object.size(merged), units = "GB")

# remove stale layers of scale.data
stale_layers <- grep("^scale\\.data\\.", Layers(integrated), value = TRUE)
for (lyr in stale_layers) {
  LayerData(integrated, assay = "RNA", layer = lyr) <- NULL
}

# same for merged if needed
stale_layers <- grep("^scale\\.data\\.", Layers(merged), value = TRUE)
for (lyr in stale_layers) {
  LayerData(merged, assay = "RNA", layer = lyr) <- NULL
}

# format(object.size(integrated), units = "GB")
# format(object.size(merged), units = "GB")

# -------------------------------------------------------------------------

# ======================================================================
# == save merged output and free up RAM ==
# ======================================================================

# save the merged object
saveRDS(merged, out_merged_object)

# free up RAM
remove(merged)
gc()

# ======================================================================
# == markers on integrated object ==
# ======================================================================

integrated.markers = lapply(resolutions, function(res) {
  # set identity to the desired clustering resolution
  Idents(integrated) = integrated[[paste0(DefaultAssay(integrated), '_snn_res.',res)]][[1]]
  markers.res = FindAllMarkers(
    integrated, 
    assay = NULL,
    features = NULL,
    logfc.threshold = 0.1, # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells
    test.use = "wilcox",
    slot = "data",
    min.pct = min_pct, # SdP: default is 0.01
    min.diff.pct = -Inf, # only test genes that show a minimum difference in the fraction of detection between the two groups.
    node = NULL, # A node to find markers for and all its children; requires BuildClusterTree to have been run previously; replaces FindAllMarkersNode
    verbose = TRUE,
    only.pos = TRUE, # SdP: default is FALSE
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
# == save integrated outputs ==
# ======================================================================

saveRDS(integrated, out_integrated_object)
write_tsv(integrated@meta.data %>% 
            rownames_to_column("barcodes"), out_meta)
write_tsv(integrated.markers, out_markers)
