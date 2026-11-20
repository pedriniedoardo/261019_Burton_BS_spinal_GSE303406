# AIM ---------------------------------------------------------------------
# on the in-house brainsphere (BS) dataset, make the same outputs as 00_test_add_metadata.R:
# proportion per sample of the expertAnno.l1 annotation and the marker panel dotplot
# (single object, single annotation: no default vs cellbender, no cytetype vs sctype)

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(cowplot)
library(scCustomize)

# custom functions --------------------------------------------------------
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

# proportion per sample of one annotation: table + boxplot, x labels carry the total cell count
prop_by_sample <- function(sobj, annotation_col, sample_col = "orig.ident", title) {
  meta <- sobj@meta.data

  # total cells per annotation (all samples) and order by overall proportion
  df_total <- meta %>%
    dplyr::count(.data[[annotation_col]], name = "n_total") %>%
    rename(annotation = 1) %>%
    mutate(prop_total = n_total / sum(n_total)) %>%
    arrange(prop_total) %>%
    mutate(label = paste0(annotation, "\n(n=", n_total,", p=", round(prop_total*100,0),"%)"))

  df_prop <- meta %>%
    group_by(.data[[sample_col]], .data[[annotation_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    rename(sample = 1, annotation = 2) %>%
    # add the zero counts for the sample x annotation combinations without cells
    complete(sample, annotation, fill = list(n = 0)) %>%
    group_by(sample) %>%
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

# output directories ------------------------------------------------------
dir_tab <- "../../out/tables/R_analysis_R45_pixi"
dir_plot <- "../../out/plots/R_analysis_R45_pixi"

# read in the data --------------------------------------------------------
sobj_bs <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_batch_2023/out/object/sobj_processed_donor.rds")

col_anno <- "expertAnno.l1"
col_sample <- "harmonized_donor2"

# filter out the doublets and the unassigned cells (donor label missing, "doublet" or "unassigned")
# check the donor labels before the filtering
table(sobj_bs@meta.data[[col_sample]], useNA = "ifany")

cells_keep <- sobj_bs@meta.data %>%
  rownames_to_column("barcode_id") %>%
  filter(!is.na(.data[[col_sample]]),
         !(.data[[col_sample]]) %in% c("doublet", "unassigned")) %>%
  pull(barcode_id)

message("cells removed: ", ncol(sobj_bs) - length(cells_keep), " of ", ncol(sobj_bs))
sobj_bs <- subset(sobj_bs, cells = cells_keep)

# check the donor labels after the filtering
table(sobj_bs@meta.data[[col_sample]], useNA = "ifany")

# confirm the annotation
DimPlot_scCustom(sobj_bs, group.by = col_anno, label = TRUE)
table(sobj_bs@meta.data[[col_sample]], sobj_bs@meta.data[[col_anno]], useNA = "ifany")

# umap --------------------------------------------------------------------
# umap only (no legend), labels with a white halo as in 00_test_add_metadata.R
p_umap <- dimplot_halo(sobj_bs, col_anno) +
  NoLegend() + ggtitle(paste("BS -", col_anno))
p_umap
ggsave(plot = p_umap,
       filename = file.path(dir_plot, "01_BS_p_umap_expertAnno_l1.pdf"),
       width = 6, height = 5)

# proportion per sample ---------------------------------------------------
res_prop <- prop_by_sample(sobj_bs, col_anno, col_sample, title = "BS")
res_prop$plot
write_csv(res_prop$table, file.path(dir_tab, "01_BS_prop_expertAnno_l1_by_donor.csv"))
ggsave(plot = res_prop$plot,
       filename = file.path(dir_plot, "01_BS_p_prop_expertAnno_l1_by_donor.pdf"),
       width = 9, height = 6)

# dotplot of the marker panel ---------------------------------------------
geneset <- list(
  IMMUNE = c("AIF1", "TYROBP", "HLA-DRA", "TREM2", "CX3CR1", "C3", "CSF1R", "CD74", "RUNX1", "C1QB", "PTPRC"),
  OLIGOLINEAGE = c("CSPG4", "OLIG1", "OLIG2", "PDGFRA", "SOX6", "PLP1", "SOX10", "BCAS1", "MBP", "MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM", "SLC1A2", "S100B"),
  NEURONS = c("GAD2", "TLE4", "CUX2", "SYP", "NEFL", "SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

# genes missing from the object would be silently dropped by DotPlot: report them
print(setdiff(unlist(geneset), rownames(sobj_bs)))

p_dotplot <- plot_panel_dotplot(sobj_bs, col_anno) + ggtitle(paste("BS -", col_anno))
p_dotplot
ggsave(plot = p_dotplot,
       filename = file.path(dir_plot, "01_BS_p_dotplot_panel_expertAnno_l1.pdf"),
       width = 16, height = 5)
