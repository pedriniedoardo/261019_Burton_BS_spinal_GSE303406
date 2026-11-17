# AIM ---------------------------------------------------------------------
# add the metadata from geo to the reanalyzed object

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(magick)
library(patchwork)
library(grid)
library(cowplot)
library(scCustomize)

# custom functions --------------------------------------------------------
# jaccard index between two sets of barcodes
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

# build the dataset crossing the two annotations per cell, and compute the jaccard score per pair
compute_jaccard_matrix <- function(meta, id_ref_col, id_query_col) {
  df_crossing <- crossing(id_ref = unique(meta[[id_ref_col]]),
                          id_query = unique(meta[[id_query_col]]))
  
  df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                                id_query = df_crossing$id_query), function(id_ref, id_query){
                                  
                                  a <- meta %>%
                                    filter(.data[[id_ref_col]] == id_ref) %>% pull(barcode)
                                  
                                  b <- meta %>%
                                    filter(.data[[id_query_col]] == id_query) %>% pull(barcode)
                                  
                                  jaccard_score <- jaccard(a, b)
                                  
                                  data.frame("id_ref" = id_ref,
                                             "id_query" = id_query,
                                             "jaccard_score" = jaccard_score)
                                }) %>%
    bind_rows()
  
  # shape it as a matrix
  df_jaccard_score %>%
    pivot_wider(names_from = id_ref, values_from = jaccard_score) %>%
    column_to_rownames("id_query")
}

# DimPlot_scCustom without its own labels, plus ggrepel labels at the median of each group with a white halo
dimplot_halo <- function(sobj, group_col, label_size = 3, bg_color = "white", bg_r = 0.15, ...) {
  p <- DimPlot_scCustom(sobj, group.by = group_col, label = FALSE, ...)

  # the first two columns of the plot data are the embedding dimensions
  dims <- colnames(p$data)[1:2]
  label_data <- p$data %>%
    mutate(label = .data[[group_col]]) %>%
    group_by(label) %>%
    summarise(x = median(.data[[dims[1]]]), y = median(.data[[dims[2]]]), .groups = "drop")

  p +
    ggrepel::geom_text_repel(data = label_data,
                             mapping = aes(x = x, y = y, label = label),
                             inherit.aes = FALSE,
                             size = label_size,
                             fontface = "bold",
                             bg.color = bg_color,
                             bg.r = bg_r,
                             max.overlaps = Inf,
                             seed = 1)
}

# full comparison for one object: umaps of both annotations + jaccard heatmap
col_auto <- "cytetype_cellOntologyTerm_RNA_snn_res.0.3"
col_authors <- "sctype_classification"

plot_annotation_comparison <- function(sobj, title) {
  # drop cells without an author annotation (barcodes not in LUT_authors)
  meta <- sobj@meta.data %>%
    rownames_to_column("barcode") %>%
    filter(!is.na(.data[[col_authors]]), !is.na(.data[[col_auto]]))
  
  # rows: automatic annotation, columns: authors' annotation
  mat_jaccard <- compute_jaccard_matrix(meta, col_authors, col_auto)
  
  ht_jaccard <- Heatmap(mat_jaccard,
                        name = "Jaccard score",
                        col = viridis::viridis(option = "turbo", n = 20),
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize = 8),
                        column_names_side = "bottom",
                        column_names_gp = gpar(fontsize = 8),
                        row_dend_reorder = FALSE,
                        column_dend_reorder = FALSE,
                        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                        column_title = "authors annotation (sctype_classification)",
                        row_title = "automatic annotation (cytetype, res 0.3)",
                        show_column_names = TRUE,
                        show_row_names = TRUE)
  
  p01 <- dimplot_halo(sobj, col_auto) +
    NoLegend() + ggtitle("cytetype (automatic)")
  
  p02 <- dimplot_halo(sobj, col_authors) +
    NoLegend() + ggtitle("sctype (authors)")
  
  # grid.grabExpr(draw(x)) turns the ComplexHeatmap into a grob that patchwork can place
  pfinal <- (p01 + p02 + grid.grabExpr(draw(ht_jaccard))) +
    plot_layout(widths = c(1, 1, 2)) +
    plot_annotation(title = title)
  
  return(pfinal)
}

# dotplot ordered by the marker group, one row per group of the chosen annotation
plot_panel_dotplot <- function(sobj, group_col, panel = geneset) {
  # source plot: only used to extract the data (clustered idents)
  p_src <- DotPlot(sobj,
                   features = panel,
                   dot.scale = 8,
                   cluster.idents = TRUE,
                   group.by = group_col)
  
  # map every feature to its marker group
  df_plot <- lapply(panel, function(x) {
    p_src$data %>% filter(features.plot %in% x)
  }) %>%
    bind_rows(.id = "marker_group") %>%
    mutate(marker_group = factor(marker_group, levels = names(panel)),
           features.plot = factor(features.plot, levels = unique(unlist(panel))))
  
  p <- ggplot(df_plot, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, col = avg.exp.scaled)) +
    scale_radius(range = c(0, 8)) +
    facet_grid(~marker_group, scales = "free", space = "free") +
    theme_cowplot() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(hjust = 1, angle = 90),
          strip.text.x = element_text(angle = 0)) +
    scale_color_gradient(low = "lightgrey", high = "blue", limits = c(-1, 2), oob = scales::squish) +
    labs(title = group_col, x = NULL, y = NULL)
  
  return(p)
}

# proportion per sample of one annotation: table + boxplot, x labels carry the total cell count
prop_by_sample <- function(sobj, annotation_col, title) {
  meta <- sobj@meta.data
  
  # total cells per annotation (all samples) and order by overall proportion
  df_total <- meta %>%
    dplyr::count(.data[[annotation_col]], name = "n_total") %>%
    rename(annotation = 1) %>%
    mutate(prop_total = n_total / sum(n_total)) %>%
    arrange(prop_total) %>%
    mutate(label = paste0(annotation, "\n(n=", n_total,", p=", round(prop_total*100,0),"%)"))
  
  df_prop <- meta %>%
    group_by(orig.ident, .data[[annotation_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    rename(annotation = 2) %>%
    # add the zero counts for the sample x annotation combinations without cells
    complete(orig.ident, annotation, fill = list(n = 0)) %>%
    group_by(orig.ident) %>%
    mutate(tot = sum(n),
           prop = n / tot) %>%
    ungroup() %>%
    left_join(df_total, by = "annotation") %>%
    mutate(label = factor(label, levels = rev(df_total$label)))
  
  p <- ggplot(df_prop, aes(x = label, y = prop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2, height = 0), shape = 1) +
    theme_cowplot() +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    labs(title = title, subtitle = annotation_col, x = NULL, y = "proportion per sample")
  
  list(table = df_prop, plot = p)
}

# output directories ------------------------------------------------------
dir_obj <- "../../out/objects/R_analysis_R45_pixi"
dir_tab <- "../../out/tables/R_analysis_R45_pixi"
dir_plot <- "../../out/plots/R_analysis_R45_pixi"

# read in teh data --------------------------------------------------------
# read in the cellbender and the default processed dataset
sobj_default <- readRDS("../../data/results/Seurat_base/object/integrated_obj_default_annotated.rds")
sobj_cellbender <- readRDS("../../data/results/Seurat_base/object/integrated_obj_cellbender_annotated.rds")

# confirm the identity of the datasets
DimPlot(sobj_default,group.by = "cytetype_cellOntologyTerm_RNA_snn_res.0.3")
DimPlot(sobj_cellbender,group.by = "cytetype_cellOntologyTerm_RNA_snn_res.0.3")

sobj_default@meta.data

# read in the metadata shared by the authors
LUT_authors <- read_csv("../../data/GSM9125398_meta_data.csv") %>%
  select(sctype_classification,cell_barcode,Sample) %>%
  mutate(barcode_id = paste(Sample,cell_barcode,sep = "_"))

meta_full_default <- sobj_default@meta.data %>%
  rownames_to_column("barcode_id") %>%
  left_join(LUT_authors,by = "barcode_id") %>%
  mutate(sctype_classification = case_when(is.na(sctype_classification)~"missing",
                                           T~sctype_classification)) %>%
  column_to_rownames("barcode_id")

meta_full_cellbender <- sobj_cellbender@meta.data %>%
  rownames_to_column("barcode_id") %>%
  left_join(LUT_authors,by = "barcode_id") %>%
  mutate(sctype_classification = case_when(is.na(sctype_classification)~"missing",
                                           T~sctype_classification)) %>%
  column_to_rownames("barcode_id")

sobj_default_full <- AddMetaData(sobj_default,metadata = meta_full_default)
sobj_cellbender_full <- AddMetaData(sobj_cellbender,metadata = meta_full_cellbender)

# confirm the addition of the metadata
DimPlot_scCustom(sobj_default_full,group.by = "cytetype_cellOntologyTerm_RNA_snn_res.0.3") | DimPlot_scCustom(sobj_default_full,group.by = "sctype_classification")
DimPlot_scCustom(sobj_cellbender_full,group.by = "cytetype_cellOntologyTerm_RNA_snn_res.0.3") | DimPlot_scCustom(sobj_cellbender_full,group.by = "sctype_classification")

# compare the automatic annotation (cytetype) vs the authors' annotation (sctype) --------
# see the implementation in the custom functions
p_default <- plot_annotation_comparison(sobj_default_full, "default processing")
ggsave(plot = p_default,filename = file.path(dir_plot,"00_p_default.pdf"),width = 15,height = 5)

p_cellbender <- plot_annotation_comparison(sobj_cellbender_full, "cellbender processing")
ggsave(plot = p_cellbender,filename = file.path(dir_plot,"00_p_cellbender.pdf"),width = 15,height = 5)

# EDA ---------------------------------------------------------------------
# estimate the prop per sample fo the OLIGO annotation

# table(sobj_default_full@meta.data$Sample,sobj_default_full@meta.data$orig.ident)

df_missing_default <- sobj_default_full@meta.data %>%
  filter(sctype_classification == "missing") %>%
  group_by(orig.ident) %>%
  summarise(n = n())
# write_csv(df_missing, file.path(dir_tab, "00_n_missing_authors_annotation_default.csv"))

df_missing_cellbender <- sobj_cellbender_full@meta.data %>%
  filter(sctype_classification == "missing") %>%
  group_by(orig.ident) %>%
  summarise(n = n())


order_id_default <- sobj_default_full@meta.data %>%
  # filter(!is.na(sctype_classification)) %>%
  group_by(sctype_classification) %>%
  summarise(n = n(),.groups = "drop") %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  arrange(prop) %>%
  pull(sctype_classification)

order_id_cellbender <- sobj_cellbender_full@meta.data %>%
  # filter(!is.na(sctype_classification)) %>%
  group_by(sctype_classification) %>%
  summarise(n = n(),.groups = "drop") %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  arrange(prop) %>%
  pull(sctype_classification)

# sobj_default_full@meta.data %>%
#   # filter(!is.na(sctype_classification)) %>%
#   group_by(orig.ident,sctype_classification) %>%
#   summarise(n = n(),.groups = "drop") %>%
#   group_by(orig.ident) %>%
#   mutate(tot = sum(n)) %>%
#   mutate(prop = n/tot) %>%
#   mutate(sctype_classification = factor(sctype_classification,levels = rev(order_id))) %>%
#   filter(sctype_classification == "missing")

list_prop_obj <- list(default = sobj_default_full, cellbender = sobj_cellbender_full)
list_prop_ann <- list(sctype = col_authors, cytetype = col_auto)

# one table + one plot per object x annotation
iwalk(list_prop_obj, function(sobj, name) {
  iwalk(list_prop_ann, function(annotation_col, annotation_name) {
    res <- prop_by_sample(sobj, annotation_col, title = name)
    print(res$plot)
    write_csv(res$table, file.path(dir_tab, paste0("00_prop_", annotation_name, "_by_sample_", name, ".csv")))
    ggsave(plot = res$plot,
           filename = file.path(dir_plot, paste0("00_p_prop_", annotation_name, "_by_sample_", name, ".pdf")),
           width = 9, height = 6)
  })
})

# dotplot of the marker panel, per annotation ------------------------------
geneset <- list(
  IMMUNE = c("AIF1", "TYROBP", "HLA-DRA", "TREM2", "CX3CR1", "C3", "CSF1R", "CD74", "RUNX1", "C1QB", "PTPRC"),
  OLIGOLINEAGE = c("CSPG4", "OLIG1", "OLIG2", "PDGFRA", "SOX6", "PLP1", "SOX10", "BCAS1", "MBP", "MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM", "SLC1A2", "S100B"),
  NEURONS = c("GAD2", "TLE4", "CUX2", "SYP", "NEFL", "SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

# compare the two annotations on the same panel, for each object
p_dotplot_default <- (plot_panel_dotplot(sobj_default_full, col_auto) / plot_panel_dotplot(sobj_default_full, col_authors)) + plot_annotation(title = "default")
p_dotplot_default
ggsave(plot = p_dotplot_default,filename = file.path(dir_plot, "00_p_dotplot_panel_default.pdf"),width = 16,height = 14)

p_dotplot_cellbender <- (plot_panel_dotplot(sobj_cellbender_full, col_auto) / plot_panel_dotplot(sobj_cellbender_full, col_authors)) + plot_annotation(title = "cellbender")
p_dotplot_cellbender
ggsave(plot = p_dotplot_cellbender,filename = file.path(dir_plot, "00_p_dotplot_panel_cellbender.pdf"),width = 16,height = 14)

# save the outputs --------------------------------------------------------
# annotated objects (with the authors' metadata)
list_sobj <- list(default = sobj_default_full,
                  cellbender = sobj_cellbender_full)

iwalk(list_sobj, function(sobj, name) {
  saveRDS(sobj, file = file.path(dir_obj, paste0("integrated_obj_", name, "_annotated_authors_meta.rds")))
})
