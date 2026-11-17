# ======================================================================
# == load libraries ==
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
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

qc_features      = snakemake@params$qc_features
figure_extension = snakemake@params$figure_extension
out_location     = snakemake@params$out_location
resolutions      = snakemake@params$resolutions
cellbender_tag   = snakemake@params$cellbender_tag
# annot_col_names  = snakemake@params$ann_cols

merge_plot_dir <- paste0(out_location, "plot/Merge/")
dir.create(paste0(out_location, "plot/UMAPs"),    showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(out_location, "plot/Barplots"), showWarnings = FALSE, recursive = TRUE)
dir.create(merge_plot_dir, recursive = TRUE, showWarnings = FALSE)

message('------------------------------------')
message("## Loading merged and integrated object...")
message('------------------------------------')
message()


pre_qc_rds        = snakemake@input$pre_qc_rds
merged_rds        = snakemake@input$merged_rds
integrated_rds    = snakemake@input$integrated_rds
samples_meta_file = snakemake@input$samples_meta
LUT_QC_file       = snakemake@input$LUT_QC

pre_qc_list    = lapply(pre_qc_rds, readRDS)
merged         = readRDS(merged_rds)
integrated     = readRDS(integrated_rds)
n_samples      = length(unique(integrated$orig.ident))
n_res          = length(resolutions)

tab_qc         = snakemake@output$tab_qc
report_done    = snakemake@output$done

message('------------------------------------')
message("## Loading samples metadata...")
message('------------------------------------')
message()

samples_meta = read.csv(samples_meta_file)

# -- restrict to the samples that survive LUT_QC's include_sample filter --
# samples_meta is hand-edited and lists every sample in samples.csv, but
# merged/integrated only contain samples with include_sample = TRUE (see
# 03_integrate.R). Matching that filter here keeps the join below exact,
# regardless of when samples_meta or LUT_QC were last edited.
df_lut = read.csv(LUT_QC_file)
included_samples = df_lut$sample_name[df_lut$include_sample]
samples_meta = samples_meta[samples_meta$sample_name %in% included_samples, ]

id_col = colnames(samples_meta)=='sample_name'
annot_col = !id_col
annot_col_names = colnames(samples_meta)[annot_col]
# rename the id colname
colnames(samples_meta)[id_col] = 'orig.ident'
# show in the log the annotation columns
message('columnn : (values)')
for( k in colnames(samples_meta)[annot_col] ) {
  message(k, ' : ', paste(samples_meta[k], collapse = ','))
}

# -------- add samples meta to the merged object ------ #

add_samples_meta = function(merged, samples_meta) {
  # move rownames into a column
  single_cell_meta = merged@meta.data %>% 
    tibble::rownames_to_column('barcode')
  # merge the two metadata (left join: keep every cell in the Seurat object,
  # drop any samples_meta rows that don't match one of its orig.ident values
  # -- avoids crashing on stray/unmatched rows regardless of their source)
  single_cell_meta =
    merge(
      single_cell_meta,
      samples_meta,
      by = 'orig.ident',
      all.x = T, sort = F) %>%
    tibble::column_to_rownames('barcode')
  # add the metadata annotations (be sure that the order is preserved)
  merged@meta.data = single_cell_meta[colnames(merged),]
  return(merged)
}

merged     = add_samples_meta(merged, samples_meta)
integrated = add_samples_meta(integrated, samples_meta)

# ======================================================================
# == set panel rows and cols ===
# ======================================================================

panel_layout <- function(n) {
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  return(c(nrow, ncol))
}

# wrap text at whitespace so long labels use several lines (vectorised)
wrap_text <- function(x, width) {
  vapply(strwrap(x, width = width, simplify = FALSE),
         paste, collapse = "\n", character(1))
}

# Save a DimPlot split by a metadata column on a roughly square grid whose canvas
# grows with the number of panels, so the plot stays readable for any number of
# samples/levels (Seurat's default is one row on a fixed canvas).
#   file_stub    output path without extension
#   exts         file extensions to write (builds the plot once for all of them)
#   raster       TRUE: points drawn as one bitmap per panel (scattermore) inside
#                the file, axes/text stay vector -> much smaller PDFs.
#                FALSE: every point is a vector shape (Seurat's geom_point)
#   pt_size      point size (NULL: 3 px when raster, Seurat's auto size otherwise)
#   legend_wrap  wrap legend labels / strip titles at this many characters, and
#                size the canvas to the widest wrapped label so the UMAP is not
#                squeezed by a long legend
save_split_dimplot <- function(file_stub, object, split_by, group_by = NULL,
                               exts = figure_extension, panel_size = 3.5,
                               max_cols = 5, raster = TRUE, title = NULL,
                               cols = NULL, legend_wrap = 30, pt_size = NULL, ...) {
  n = length(unique(stats::na.omit(object@meta.data[[split_by]])))
  if (n < 1) {
    message("skipping ", basename(file_stub), ": '", split_by, "' has no non-NA values")
    return(invisible(NULL))
  }
  ncols = min(ceiling(sqrt(n)), max_cols)
  nrows = ceiling(n / ncols)

  # Seurat's raster default is a 1-pixel point (too faint on a ~700 px panel);
  # for vector output keep Seurat's auto size
  if (is.null(pt_size) && raster) pt_size = 3
  p = DimPlot(object, group.by = group_by, split.by = split_by, ncol = ncols,
              cols = cols, raster = raster,
              pt.size = pt_size,
              raster.dpi = rep(round(panel_size * 200), 2),  # bitmap pixels per panel
              ...)

  # legend labels: wrap long ones (re-declare the colour scale, keeping `cols`)
  lg_lines = 0; lg_chars = 0
  if (!is.null(group_by) && group_by %in% colnames(object@meta.data)) {
    lv = as.character(unique(stats::na.omit(object@meta.data[[group_by]])))
    wrapped  = wrap_text(lv, legend_wrap)
    lg_lines = sum(lengths(strsplit(wrapped, "\n", fixed = TRUE)))
    lg_chars = max(nchar(unlist(strsplit(wrapped, "\n", fixed = TRUE))), 0)
    lab_fun  = function(x) wrap_text(x, legend_wrap)
    p = suppressMessages(
      p + if (is.null(cols)) scale_colour_discrete(labels = lab_fun)
          else scale_colour_manual(values = cols, labels = lab_fun))
  }

  # long panel (split level) titles: wrap them too
  if (!is.null(p$facet$params$labeller))
    p$facet$params$labeller = ggplot2::label_wrap_gen(width = legend_wrap, multi_line = TRUE)

  # title: left-aligned on the whole plot (a long title is no longer clipped)
  p = p + labs(title = if (is.null(title)) group_by else title) +
    theme(strip.text         = element_text(size = 8),
          legend.text        = element_text(size = 7),
          legend.key.size    = grid::unit(0.8, "lines"),
          plot.title         = element_text(size = 12, face = "bold", hjust = 0),
          plot.title.position = "plot")

  # fix the legend layout ourselves (ggplot would wrap it into columns on its own and
  # squeeze the panels): as many columns as needed so every column fits the plot height
  lg_line_h = 0.2                                # inches per legend line
  max_rows  = max(1, floor((nrows * panel_size + 0.7) / lg_line_h))
  lg_ncol   = max(1, ceiling(lg_lines / max_rows))
  if (lg_lines > 0)
    p = p + guides(colour = guide_legend(ncol = lg_ncol, byrow = FALSE,
                                         override.aes = list(size = 2.5)))

  legend_w = lg_ncol * (0.35 + 0.06 * lg_chars) + 0.3   # inches: key + label per column
  width  = ncols * panel_size + legend_w
  height = nrows * panel_size + 0.7
  for (ext in exts) {
    ggsave(paste0(file_stub, ".", ext), p,
           width = width, height = height, limitsize = FALSE)
  }
}

# Save a (non-split) DimPlot with repelled labels and a compact legend whose
# canvas follows the number/length of the legend entries.
#   label  draw the group names on the clusters (repel = TRUE avoids overlaps)
save_dimplot <- function(file_stub, object, group_by, exts = figure_extension,
                         label = FALSE, title = NULL, panel_size = 6,
                         legend_wrap = 30, legend_text_size = 7, ...) {
  p = DimPlot(object, group.by = group_by, label = label, repel = TRUE,
              label.size = 3, ...)
  p = p + labs(title = if (is.null(title)) group_by else title) +
    theme(legend.text         = element_text(size = legend_text_size),
          legend.key.size     = grid::unit(0.8, "lines"),
          plot.title          = element_text(size = 12, face = "bold", hjust = 0),
          plot.title.position = "plot")

  lv = as.character(unique(stats::na.omit(object@meta.data[[group_by]])))
  n_lv = length(lv)
  if (!is.numeric(object@meta.data[[group_by]])) {
    p = suppressMessages(p + scale_colour_discrete(labels = function(x) wrap_text(x, legend_wrap)))
    lines = unlist(strsplit(wrap_text(lv, legend_wrap), "\n", fixed = TRUE))
    chars = max(nchar(lines), 0)
    # legend columns: ~30 entries per column, as ggplot would wrap them
    n_col = ceiling(length(lines) / 30)
    p = p + guides(colour = guide_legend(ncol = n_col, override.aes = list(size = 2.5)))
    legend_w = n_col * (0.5 + 0.06 * chars)
  } else legend_w = 1

  width  = panel_size + 0.6 + legend_w
  height = panel_size + 0.7
  for (ext in exts) {
    ggsave(paste0(file_stub, ".", ext), p,
           width = width, height = height, limitsize = FALSE)
  }
}

n_rows_clusters = panel_layout(n_res)[1]
n_cols_clusters = panel_layout(n_res)[2]

n_rows_samples  = panel_layout(n_samples)[1]
n_cols_samples  = panel_layout(n_samples)[2]

message('------------------------------------')
message("## merge preQC objects...")
message('------------------------------------')
message()

# to avoid the warning when merging:
pre_qc_seurat_list = lapply(pre_qc_list, function(obj) {
  RenameCells(obj, add.cell.id = Project(obj)) 
  })
# merge actually
pre_qc_merged = merge(pre_qc_seurat_list[[1]], pre_qc_seurat_list[-1])

message('------------------------------------')
message("## QC violin plots...")
message('------------------------------------')
message()

# remove qc features not computed in the preQC objects
pre_qc_features = intersect(qc_features, colnames(pre_qc_merged@meta.data))

# plot actually: canvas grows with the number of samples (3 panels per row),
# sample names are rotated and shrunk, no per-sample legend
n_vln      = length(pre_qc_features)
vln_ncol   = min(3, n_vln)
vln_nrow   = ceiling(n_vln / vln_ncol)
vln_panel_w = max(4, 0.3 * n_samples + 1.5)
p_vln_qc = VlnPlot(pre_qc_merged, log = T,
                   layer = 'counts',
                   pt.size = 0,
                   features = pre_qc_features,
                   group.by = 'orig.ident',
                   ncol = vln_ncol) &
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none")
for(ext in figure_extension) {
  ggsave(
    paste0(out_location, "plot/vln_", cellbender_tag, ".", ext),
    p_vln_qc,
    height = vln_nrow * (3.5 + 0.055 * max(nchar(unique(as.character(pre_qc_merged$orig.ident))))),
    width  = vln_ncol * vln_panel_w,
    limitsize = FALSE)
}

message('------------------------------------')
message("## QC stats...")
message('------------------------------------')
message()

message('1')

pre_filtering_stats = 
  pre_qc_merged@meta.data %>% 
  group_by(orig.ident) %>% 
  summarise(
    n_cells = n(),
    mean_UMI = round(mean(nCount_RNA)),
    mean_genes = round(mean(nFeature_RNA))
  )

message('2')

post_filtering_stats = 
  merged@meta.data %>% 
  group_by(orig.ident) %>% 
  summarise(
    n_cells = n(),
    mean_UMI = round(mean(nCount_RNA)),
    mean_genes = round(mean(nFeature_RNA))
  )

message('3')

QC_stats = 
  merge(
    pre_filtering_stats, 
    post_filtering_stats, 
    by = 'orig.ident', 
    suffixes = c('(preQC)', '(postQC)')
  )

message('4')

write_csv(QC_stats, file = tab_qc)

message('------------------------------------')
message("## UAMP of QC features...")
message('------------------------------------')
message()

for(feat in qc_features) {
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "plot/UMAPs/dim_", feat, "_", cellbender_tag, ".", ext),
      FeaturePlot(integrated, features = feat), 
      height = 6, width = 8
      )
  }
}

message('------------------------------------')
message("## UMAP of samples...")
message('------------------------------------')
message()

# merged and integrated

for(object_type in c('integrated','merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
                  )
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "plot/UMAPs/dim_", object_type, "_ident_", cellbender_tag, ".", ext),
      DimPlot(object, group.by = 'orig.ident'),
      height = 6, width = 8
    )
  }
}

message('------------------------------------')
message("## UMAP of clusters...")
message('------------------------------------')
message()

# merged and integrated

for(object_type in c('integrated','merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
  )
  for(ext in figure_extension) {
    for( res in resolutions ) {
      ggsave(
        paste0(out_location, "plot/UMAPs/dim_", object_type, "_cluster_RES", 
               res, "_", cellbender_tag, ".", ext),
        DimPlot(object, group.by = paste0('RNA_snn_res.', res)),
        height = 6, width = 8)
    }
  }
}

# split by annotations of the samples

for(ann in annot_col_names) {
  for(object_type in c('integrated','merged')) {
    object = switch(object_type, 
                    'integrated' = integrated, 
                    'merged'     = merged
    )
    for( res in resolutions ) {
      save_split_dimplot(
        paste0(out_location, "plot/UMAPs/dim_", object_type, "_cluster_RES",
               res, "_by", ann, "_", cellbender_tag),
        object,
        split_by = ann,
        group_by = paste0('RNA_snn_res.', res))
    }
  }
}

message('------------------------------------')
message("## UMAP of celltype annotations")
message('------------------------------------')
message()

# merged and integrated

# retrieve the annotation metadata: SingleR's CELLTYPE_ prefix, plus
# CyteTypeR's annotation/cellOntologyTerm/cellState columns (skip
# cellOntologyTermID and granularAnnotation -- not interesting to plot).
# CyteTypeR has no numeric confidence/quality score, unlike SingleR's
# delta.next (CELLTYPEQ_ below), so it only ever contributes to this
# categorical DimPlot block.
celltype_meta = grep(
  'CELLTYPE_|cytetype_annotation_|cytetype_cellOntologyTerm(?!ID)|cytetype_cellState_',
  colnames(integrated@meta.data),
  value = T, perl = T
)

for(object_type in c('integrated')) { #,'merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
  )
  for(ct_meta in celltype_meta) {
    save_dimplot(
      paste0(out_location, "plot/UMAPs/dim_", object_type, "_", ct_meta,
             "_", cellbender_tag),
      object, group_by = ct_meta, label = TRUE)
  }
}

# split by annotations of the samples

for(ann in annot_col_names) {
  for(object_type in c('integrated')) { #,'merged')) {
    object = switch(object_type, 
                    'integrated' = integrated, 
                    'merged'     = merged
    )
    for(ct_meta in celltype_meta) {
      save_split_dimplot(
        paste0(out_location, "plot/UMAPs/dim_", object_type, "_",ct_meta,
               "_by", ann, "_", cellbender_tag),
        object,
        split_by = ann,
        group_by = ct_meta)
    }
  }
}

message('------------------------------------')
message("## UMAP of celltype annotations qualities")
message('------------------------------------')
message()

# merged and integrated

# retrive the annotation metadata thanks to the preifx CELLTYPE_
celltypeq_meta = grep('CELLTYPEQ_', 
                     colnames(integrated@meta.data), 
                     value = T)

for(object_type in c('integrated')) { #,'merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
  )
  for(ext in figure_extension) {
    for(ctq_meta in celltypeq_meta) {
      ggsave(
        paste0(out_location, "plot/UMAPs/dim_", object_type, "_", ctq_meta,  
               "_", cellbender_tag, ".", ext),
        FeaturePlot(object, features = ctq_meta),
        height = 6, width = 8)
    }
  }
}


message('------------------------------------')
message("## CellCycle dim plots...")
message('------------------------------------')
message()

# merged and integrated

for(object_type in c('integrated','merged')) {
  object = switch(object_type, 
                  'integrated' = integrated, 
                  'merged'     = merged
  )
  for(ext in figure_extension) {
      ggsave(
        paste0(out_location, "plot/UMAPs/dim_", object_type, "_CellCycle_", 
               cellbender_tag, ".", ext),
        DimPlot(object, group.by = "Phase"),
        height = 6, width = 8)
  }
}

for(object_type in c('integrated','merged')) {
  object = switch(object_type,
                  'integrated' = integrated,
                  'merged'     = merged
  )
  save_split_dimplot(
    paste0(out_location, "plot/UMAPs/dim_", object_type, "_CellCycle_split_by_sample",
           cellbender_tag),
    object,
    split_by = "orig.ident",
    group_by = "Phase")
}


message('------------------------------------')
message("## Moving merged plots to Merge folder...")
message('------------------------------------')

merged_files <- list.files(
  path       = paste0(out_location, "plot/UMAPs/"),
  pattern    = paste0("_merged_.*\\.(", paste(figure_extension, collapse = "|"), ")$"),
  full.names = TRUE
)

file.rename(from = merged_files, to = file.path(merge_plot_dir, basename(merged_files)))


message('------------------------------------')
message("## Additional report images...")
message('------------------------------------')
message()

# ---- cluster colour palette (Seurat default hue_pal) --------------------

cluster_colors = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# ---- load markers -------------------------------------------------------
# integrated_markers_default.tsv: columns cluster, gene, avg_log2FC, resolution
all_markers_full = read.table(snakemake@input$markers_tsv,
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- Dotplot + top5Heatmap (one per resolution) -----------------
for(res in resolutions) {
  
  res_col = paste0("RNA_snn_res.", res)
  Idents(integrated) = integrated@meta.data[[res_col]]
  
  # filter markers for this resolution
  all_markers = all_markers_full %>%
    filter(resolution == res) %>%
    mutate(cluster = as.character(cluster))
  
  top5 = all_markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    ungroup()
  top5_genes = unique(top5$gene)
  
  n_cl = length(unique(integrated@meta.data[[res_col]]))
  
  message(paste0("## Dotplot res=", res, "..."))
  
  p_dot = DotPlot(integrated,
                  features = top5_genes,
                  group.by = res_col,
                  scale    = TRUE) +
    scale_color_gradient2(
      low = "white", mid = "#e8a0b4", high = "#8b0057",
      midpoint = 0, name = "relative\nexpression") +
    scale_size_continuous(range = c(0.5, 6), name = "percent\nexpression") +
    labs(x = NULL, y = NULL,
         title = "Top 5 specific markers for each cluster") +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          axis.text.y  = element_text(size = 9),
          panel.grid   = element_line(color = "grey92"),
          plot.title   = element_text(hjust = 0.5))
  
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "plot/Dotplot_RNA_snn_res.", res, ".", ext),
      p_dot,
      height = max(6,  n_cl * 0.5 + 2),
      width  = max(14, length(top5_genes) * 0.28 + 4),
      limitsize = FALSE)
  }
  
  message(paste0("## top5Heatmap res=", res, "..."))
  
  p_heatmap = DoHeatmap(
    integrated,
    features = top5_genes,
    group.by = res_col,
    slot     = "scale.data",
    size     = 3,
    angle    = 0) +
    scale_fill_gradientn(
      colors   = c("#9400d3", "#000000", "#ffff00"),
      name     = "Expression",
      na.value = "white") +
    theme(axis.text.y  = element_text(size = 7),
          legend.title = element_text(size = 9))
  
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "plot/top5Heatmap_RNA_snn_res.", res, ".", ext),
      p_heatmap,
      height = max(10, length(top5_genes) * 0.22 + 3),
      width  = 20,
      limitsize = FALSE)
  }
}

# ---- Vln_per_cluster (one file per resolution) ----------------------
message("## Vln_per_cluster...")

vln_features = intersect(
  c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.globin"),
  colnames(integrated@meta.data))

for(res in resolutions) {
  res_col = paste0("RNA_snn_res.", res)
  n_cl    = length(unique(integrated@meta.data[[res_col]]))
  
  p_vln = VlnPlot(
    integrated,
    features = vln_features,
    group.by = res_col,
    pt.size  = 0,
    log      = FALSE,
    combine  = TRUE) &
    theme(axis.text.x  = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          plot.title   = element_text(size = 10, face = "bold"))
  
  for(ext in figure_extension) {
    ggsave(
      paste0(out_location, "plot/Vln_per_cluster_RNA_snn_res.", res, ".", ext),
      p_vln,
      height = 4 * ceiling(length(vln_features) / 2),
      width  = max(12, n_cl * 0.9),
      limitsize = FALSE)
  }
}


message("## CC_bySample...")

cc_phase_colors = c("G1" = "#f08080", "G2M" = "#4caf50", "S" = "#6baed6")

cc_df = integrated@meta.data %>%
  select(orig.ident, Phase) %>%
  group_by(orig.ident, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  mutate(Phase = factor(Phase, levels = c("S", "G2M", "G1")))

p_cc = ggplot(cc_df, aes(x = orig.ident, y = freq, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(values = cc_phase_colors, name = "Cell cycle phase") +
  labs(x = NULL, y = "Cell fraction") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

n_samples_w = max(8, length(unique(integrated$orig.ident)) * 1.4)

for(ext in figure_extension) {
  ggsave(paste0(out_location, "plot/Barplot_CC_bySample.", ext),
         p_cc, height = 7, width = n_samples_w)
}

# ---- per-resolution UMAPs and Barplots ----------------------------------
cond_col = if (length(annot_col_names) > 0) annot_col_names[1] else "orig.ident"

for(res in resolutions) {

  res_col   = paste0("RNA_snn_res.", res)
  cl_levels = as.character(sort(as.numeric(
    unique(as.character(integrated@meta.data[[res_col]])))))
  n_cl      = length(cl_levels)
  cl_pal    = setNames(cluster_colors(n_cl), cl_levels)

  integrated@meta.data[[res_col]] =
    factor(integrated@meta.data[[res_col]], levels = cl_levels)

  # UMAP split by condition
  message(paste0("## UMAP_byCond res=", res, "..."))
  n_conds = length(unique(integrated@meta.data[[cond_col]]))

  save_split_dimplot(
    paste0(out_location, "plot/UMAPs/UMAP_byCond_RNA_snn_res.", res),
    integrated,
    split_by = cond_col,
    group_by = res_col,
    exts     = "png",
    title    = paste0("RNA_snn_res.", res),
    label    = FALSE,
    cols     = cl_pal)

  # Barplot by sample
  message(paste0("## Barplot_bySample res=", res, "..."))

  bar_sample = integrated@meta.data %>%
    select(orig.ident, cluster = !!res_col) %>%
    group_by(orig.ident, cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(cluster = factor(cluster, levels = cl_levels))

  p_bar_sample = ggplot(bar_sample,
                        aes(x = orig.ident, y = freq, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = cl_pal, name = "Cluster") +
    labs(x = NULL, y = "Cell fraction") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))
  
  ggsave(
    paste0(out_location, "plot/Barplots/Barplot_bySample_RNA_snn_res.", res, ".png"),
    p_bar_sample,
    height = 7,
    width  = max(8, length(unique(integrated$orig.ident)) * 1.4))

  # Barplot by condition new
  message(paste0("## Barplot_byCond res=", res, "..."))

  bar_cond = integrated@meta.data %>%
    select(condition = !!cond_col, cluster = !!res_col) %>%
    group_by(condition, cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(condition) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(cluster = factor(cluster, levels = cl_levels))

  p_bar_cond = ggplot(bar_cond,
                      aes(x = condition, y = freq, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = cl_pal, name = "Cluster") +
    labs(x = NULL, y = "Cell fraction") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

  ggsave(
    paste0(out_location, "plot/Barplots/Barplot_byCond_RNA_snn_res.", res, ".png"),
    p_bar_cond,
    height = 7,
    width  = max(6, n_conds * 2.5))
}


message('------------------------------------')
message("## Report done...")
message('------------------------------------')
message()

write(NULL, file = report_done)
