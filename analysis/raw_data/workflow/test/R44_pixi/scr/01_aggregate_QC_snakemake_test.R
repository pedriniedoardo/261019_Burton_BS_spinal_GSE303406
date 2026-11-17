# AIM ---------------------------------------------------------------------
# make the variables loaded from the config file

# ======================================================================
# == load libraries ==
# ======================================================================

library(tidyverse)
library(scales)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# define the meta tables
# meta_tables = expand(rules.runSamplePreprocessing.output.meta,sample_name=SAMPLES.keys())
meta_tables <- c("/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test/results/Seurat_base/table/CD56_Neg_PAO12790_default_meta_preQC.tsv",
                 "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test/results/Seurat_base/table/CD56_Plus_PAO12790_default_meta_preQC.tsv",
                 "/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test/results/Seurat_base/table/Unsorted_PAO12790_default_meta_preQC.tsv")
# meta_tables <- snakemake@input$meta_tables
message("input metadata: ", paste(meta_tables,collapse = " "))

# pull all the sample names
# NOTE: Snakemake guarantees these are in the same order as meta_tables
sample_names <- c("CD56_Neg_PAO12790","CD56_Plus_PAO12790","Unsorted_PAO12790")
# sample_names <- snakemake@params$sample_names
message("sample wildcards: ", paste(sample_names,collapse = " "))

feature_filter <- c("nCount_RNA","nFeature_RNA","percent.ribo","percent.mt","percent.globin")
message("feature for filtering: ", paste(feature_filter,collapse = " "))

# define the output
# out_id_object <-  "results/Seurat/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
# out_id_meta <- "results/Seurat/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.rds"
# out_plot_mito <- "out/plot/01_test_out_histo_mito.pdf"
# out_plot_ribo <- "out/plot/01_test_out_histo_ribo.pdf"
# out_plot_hemo <- "out/plot/01_test_out_histo_hemo.pdf"
# out_plot_count <- "out/plot/01_test_out_histo_count.pdf"
# out_plot_feature <- "out/plot/01_test_out_histo_feature.pdf"
# out_LUT_QC <- "out/table/01_test_out_LUT_QC.csv"
# out_meta_total <- "out/table/01_test_out_meta_total.tsv"
# out_plot_mito <- snakemake@output$plot_mito
# out_plot_feature <- snakemake@output$plot_feature
# out_LUT_QC <- snakemake@output$LUT_QC
# out_meta_total <- snakemake@output$meta_total

# message("output plot mito: ", out_plot_mito)
# message("output plot feature: ", out_plot_feature)
# message("output table meta total before QC: ", out_LUT_QC)
# message("output table LUT QC: ", out_LUT_QC)

# make sure meta_tables and sample_names are matched
# file_path <- meta_tables[1]
# name_str <- sample_names[1]
check_alignment <- pmap(list(meta_tables, sample_names), function(file_path, name_str) {
  # 'fixed()' ensures we don't treat the name as a regex pattern
  str_detect(file_path, fixed(name_str))
}) %>%
  unlist()

# test the check
if (!all(check_alignment)) {
  mismatch_idx <- which(!check_alignment)
  stop(sprintf(
    "Order mismatch detected at indices: %s.\nFile: %s\nTarget Name: %s",
    paste(mismatch_idx, collapse = ", "),
    meta_tables[mismatch_idx[1]],
    sample_names[mismatch_idx[1]]
  ))
}
message("SUCCESS: Sample names and file paths are correctly aligned.")

# ======================================================================
# == load input ==
# ======================================================================

# define the basename
file <- basename(meta_tables)
message("input file names: ", file)

# extract the sample name keeping the tag
# sample_name <- file %>%
#   str_remove(pattern = "_meta_preQC.tsv")

# # extract the filename from the table to keep the matching
# sample_name <- file %>%
#   str_extract(pattern = paste(sample_names,collapse = "|"))
# message("input sample names: ", sample_name)

meta_total <- pmap(list(meta_tables,sample_names),function(x,x_name){
  read_tsv(x) %>%
    mutate(dataset = x_name)
}) %>%
  bind_rows()

# define some plotting parameters
# define number of rows in the panel
test <- length(sample_names)

message("number of samples processed: ", test)

# test <- 10
panel_row <- test %>% sqrt() %>% round(digits = 0)
panel_col <- (test / panel_row) %>% ceiling()

# ======================================================================
# == run the standard processing ==
# ======================================================================

# print head of the full metadata and its dimensions
head(meta_total)
dim(meta_total)

# automatic "nice" breaks for the log1p-scaled axes, derived from each panel's own data range instead of a metric-specific hardcoded vector.
# 

range

breaks_test <- function(x){
  seq(from = range(log1p(x))[1], to =range(log1p(x))[2],length.out = 9) %>%
    expm1(.) %>%
    round(digits = 0)
}
# x <- meta_total$percent.mt
# breaks_test <- seq(from = range(log1p(x))[1], to =range(log1p(x))[2],length.out = 9) %>%
#     expm1(.) %>%
#     round(digits = 0)

x <- meta_total$nCount_RNA
x <- meta_total$percent.mt
x <- meta_total$percent.globin

# make the breaks nicer
breaks_log1p <- function(x) {
  hi <- max(c(x, 0), na.rm = TRUE)
  low <- min(c(x), na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) return(0)
  unique(c(0, scales::breaks_log(n = 10)(c((low+1), hi))))
}

# plot mito prop
p_mito <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.ribo,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  # geom_vline(xintercept = c(10),col="red",linetype="dashed") +
  # annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

p_ribo <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.ribo,percent.mt,percent.globin)) %>%
  filter(var == "percent.ribo") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  # geom_vline(xintercept = c(10),col="red",linetype="dashed") +
  # annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

p_hemo <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.ribo,percent.mt,percent.globin)) %>%
  filter(var == "percent.globin") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  # geom_vline(xintercept = c(10),col="red",linetype="dashed") +
  # annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# plot counts
p_count <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nCount_RNA") %>%
  ggplot(aes(x=value)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  # geom_vline(xintercept = c(500,60000),col="red",linetype="dashed") +
  # annotate("rect", xmin=500, xmax=60000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# plot features
p_feature <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(dataset~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = breaks_log1p) +
  # geom_vline(xintercept = c(1000,6000),col="red",linetype="dashed") +
  # annotate("rect", xmin=1000, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# define the table to fill the QC parameters
# after the lates FG we decided to define some all inclusive aprameters
# LUT_df <- data.frame(sample_name = sample_name,
#                      featureLow_thr = NA,
#                      featureHigh_thr = NA,
#                      mito_thr = NA)

LUT_df <- data.frame(sample_name = sample_names,
                     nCount_RNA_high = Inf,
                     nCount_RNA_low = 0,
                     nFeature_RNA_high = Inf,
                     nFeature_RNA_low = 0,
                     percent.ribo_high = Inf,
                     percent.ribo_low = 0,
                     percent.mt_high = Inf,
                     percent.mt_low = 0,
                     percent.globin_high = Inf,
                     percent.globin_low = 0,
                     remove_doublets = F
                     )

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

# save tables
write_tsv(meta_total,out_meta_total)
write_csv(LUT_df,out_LUT_QC)

