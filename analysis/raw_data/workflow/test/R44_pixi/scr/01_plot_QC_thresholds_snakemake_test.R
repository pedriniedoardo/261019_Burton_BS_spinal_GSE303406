# ======================================================================
# == load libraries ==
# ======================================================================

library(tidyverse)
library(scales)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

meta_total_path <- "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test/results/Seurat_base/table/meta_total_beforeQC_V5_default.tsv"
# meta_total_path <- snakemake@input$meta_total
message("input meta total: ", meta_total_path)

# read in the LUT_QC
# LUT_QC_path <- snakemake@input$LUT_QC
LUT_QC_path <- "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test/results/Seurat_base/table/LUT_QC_default.csv"

# LUT_QC_path <- snakemake@input$LUT_QC
message("input LUT QC: ", LUT_QC_path)

# out_plot_mito <- snakemake@output$plot_mito
# out_plot_ribo <- snakemake@output$plot_ribo
# out_plot_hemo <- snakemake@output$plot_hemo
# out_plot_count <- snakemake@output$plot_count
# out_plot_feature <- snakemake@output$plot_feature
# 
# message("output plot mito: ", out_plot_mito)
# message("output plot feature: ", out_plot_feature)

# ======================================================================
# == load input ==
# ======================================================================

meta_total <- read_tsv(meta_total_path)
LUT_QC <- read_csv(LUT_QC_path)

# define some plotting parameters
# define number of rows in the panel
sample_names <- sort(unique(meta_total$dataset))
test <- length(sample_names)

message("number of samples processed: ", test)

panel_row <- test %>% sqrt() %>% round(digits = 0)
panel_col <- (test / panel_row) %>% ceiling()

# ======================================================================
# == run the standard processing ==
# ======================================================================

# print head of the full metadata and its dimensions
head(meta_total)
dim(meta_total)

# automatic breaks for the log1p-scaled axes, derived from each panel's own data range instead of a metric-specific hardcoded vector.
breaks_log1p <- function(x) {
  hi <- max(c(x, 0), na.rm = TRUE)
  low <- min(c(x), na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) return(0)
  unique(c(0, scales::breaks_log(n = 10)(c((low+1), hi))))
}

# reshape the (user-edited) LUT into one row per sample x metric, carrying that metric's low/high threshold
threshold_long <- LUT_QC %>%
  select(sample_name, ends_with("_low"), ends_with("_high")) %>%
  pivot_longer(
    cols = -sample_name,
    names_to = c("var", ".value"),
    names_pattern = "^(.*)_(low|high)$"
  ) %>%
  rename(dataset = sample_name)

# per-metric, per-sample shaded "kept" range (low to high)
# metric <- "percent.mt"
get_threshold_rect <- function(metric) {
  threshold_long %>% filter(var == metric)
}

# per-metric, per-sample dashed lines at the low/high bounds (dropping +-Inf defaults)
get_threshold_lines <- function(metric) {
  threshold_long %>%
    filter(var == metric) %>%
    pivot_longer(cols = c(low, high), names_to = "bound", values_to = "thr") %>%
    filter(is.finite(thr))
}

# plot mito prop
p_mito <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.ribo,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  geom_vline(data = get_threshold_lines("percent.mt"), aes(xintercept = thr), col = "red", linetype = "dashed") +
  geom_rect(data = get_threshold_rect("percent.mt"), aes(xmin = low, xmax = high), ymin = 0, ymax = Inf, inherit.aes = FALSE, alpha = 0.1, fill = "red") +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# plot ribo prop
p_ribo <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.ribo,percent.mt,percent.globin)) %>%
  filter(var == "percent.ribo") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  geom_vline(data = get_threshold_lines("percent.ribo"), aes(xintercept = thr), col = "red", linetype = "dashed") +
  geom_rect(data = get_threshold_rect("percent.ribo"), aes(xmin = low, xmax = high), ymin = 0, ymax = Inf, inherit.aes = FALSE, alpha = 0.1, fill = "red") +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# plot hemoglobin prop
p_hemo <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.ribo,percent.mt,percent.globin)) %>%
  filter(var == "percent.globin") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  geom_vline(data = get_threshold_lines("percent.globin"), aes(xintercept = thr), col = "red", linetype = "dashed") +
  geom_rect(data = get_threshold_rect("percent.globin"), aes(xmin = low, xmax = high), ymin = 0, ymax = Inf, inherit.aes = FALSE, alpha = 0.1, fill = "red") +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# plot counts
p_count <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nCount_RNA") %>%
  ggplot(aes(x=value)) +
  geom_histogram(binwidth = 0.05) +
  geom_vline(data = get_threshold_lines("nCount_RNA"), aes(xintercept = thr), col = "red", linetype = "dashed") +
  geom_rect(data = get_threshold_rect("nCount_RNA"), aes(xmin = low, xmax = high), ymin = 0, ymax = Inf, inherit.aes = FALSE, alpha = 0.1, fill = "red") +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# plot features
p_feature <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) +
  geom_histogram(binwidth = 0.05) +
  geom_vline(data = get_threshold_lines("nFeature_RNA"), aes(xintercept = thr), col = "red", linetype = "dashed") +
  geom_rect(data = get_threshold_rect("nFeature_RNA"), aes(xmin = low, xmax = high), ymin = 0, ymax = Inf, inherit.aes = FALSE, alpha = 0.1, fill = "red") +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# ======================================================================
# == save output ==
# ======================================================================

# save plot histo mito prop reads
ggsave(plot = p_mito,
       out_plot_mito,
       width = 4 * panel_col,
       height = 4 * panel_row)

# save plot histo ribo prop reads
ggsave(plot = p_ribo,
       out_plot_ribo,
       width = 4 * panel_col,
       height = 4 * panel_row)

# save plot histo globin prop reads
ggsave(plot = p_hemo,
       out_plot_hemo,
       width = 4 * panel_col,
       height = 4 * panel_row)

# save plot histo count reads
ggsave(plot = p_count,
       out_plot_count,
       width = 4 * panel_col,
       height = 4 * panel_row)

# save plot histo feature reads
ggsave(plot = p_feature,
       out_plot_feature,
       width = 4 * panel_col,
       height = 4 * panel_row)
