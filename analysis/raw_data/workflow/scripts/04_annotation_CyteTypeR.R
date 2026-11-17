# ======================================================================
# == load libraries ==
# ======================================================================

# add the installtion step of the GitHub package.
# all the dependencies have been tracked in the yaml file

# Bootstrap CyteTypeR
if (!requireNamespace("CyteTypeR", quietly = TRUE)) {
  message("Installing CyteTypeR from GitHub...")
  # fix the version
  remotes::install_github("NygenAnalytics/CyteTypeR", ref = "v0.5.1", upgrade = "never")
  
  # lates version
  # remotes::install_github("NygenAnalytics/CyteTypeR", upgrade = "never")
}

library(CyteTypeR)
library(tidyverse)
library(patchwork)
library(Matrix)
library(Seurat)
library(presto)
library(duckdb)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# input:
rds <- snakemake@input$rds
markers <- snakemake@input$markers
metadata <- snakemake@input$metadata
cellbender_tag <- snakemake@params$cellbender_tag

# rds <- "../../../results/Seurat/object/integrated_obj_cellbender.rds"
# markers <- "../../../results/Seurat/table/integrated_markers_cellbender.tsv"
# metadata <- "../../../config/meta_annotation.yaml"

message("input rds: ", rds)
message("input metadata: ", metadata)
message("cellbender tag: ", cellbender_tag)

# output:
out_id_meta <- snakemake@output$annotation_meta
out_id_varsh5 <- snakemake@output$vars_h5_path
out_id_duckdb <- snakemake@output$obs_duckdb_path
out_id_filename <- snakemake@output$query_filename
out_id_jobinfo <- snakemake@output$job_id_file

# out_id_meta <- "out/table/integrated_cytetypeR_meta_cellbender.csv"
# out_id_varsh5 <- "out/object/vars.h5"
# out_id_duckdb <- "out/object/obs.duckdb"
# out_id_filename <- "out/object/query.json"
# out_id_jobinfo <- "out/table/cytetype_job_id_cellbender.json"

message("output CyteTypeR annotation table: ", out_id_meta)
message("output CyteTypeR vars.h5: ", out_id_varsh5)
message("output CyteTypeR obs.duckdb: ", out_id_duckdb)
message("output CyteTypeR query.json: ", out_id_filename)
message("output CyteTypeR job_id file: ", out_id_jobinfo)

# params
# token <- snakemake@params$token
token <- Sys.getenv("CYTETYPER_TOKEN")
# this should be handled in the snakemake rule
if (token == "") {
  stop("CYTETYPER_TOKEN environment variable is not set or is empty. Please export it before running Snakemake.")
}

cov_markers <- snakemake@params$cov_markers

# token <- ""
# cov_markers <- "RNA_snn_res.0.3"

# ======================================================================
# == load input ==
# ======================================================================

# read in the full object
scobj <- readRDS(rds)

# read in the makrers
df_markers <- read_tsv(markers)

# read in the metadata
list_metadata <- yaml::read_yaml(metadata)

# stop early if the template placeholders were not replaced
if (any(grepl("^CHANGE_ME", unlist(list_metadata)))) {
  stop(metadata, " still contains CHANGE_ME placeholders; fill it in before running annotation.")
}

# modify the title to incorporate the cellbender tag
list_metadata$title <- paste(list_metadata$title,cellbender_tag,sep = "_")

# ======================================================================
# == run the standard processing ==
# ======================================================================

# reshape the marker output to be accepted
# This one assumes that the order of the genes is maintained from the output of Seurat's FindAllMarkers() function during the integration step
df_markers_fix <- df_markers %>%
  # mutate(cluster = as.numeric(cluster)) %>%
  filter(resolution == str_remove(cov_markers,pattern = "RNA_snn_res.")) %>%
  select(cluster,gene,avg_log2FC)

# CyteTypeR requires (1) the same set of cluster labels in the object and the marker table, and (2) >= 5 markers per cluster. FindAllMarkers drops tiny clusters (min.cells.group), so keep only the clusters satisfying both and run CyteTypeR on that subset; the dropped cells are labelled again at the end.
min_markers <- 5
obj_clusters <- as.character(unique(scobj@meta.data[[cov_markers]]))
marker_counts <- table(as.character(df_markers_fix$cluster))
keep_clusters <- intersect(obj_clusters, names(marker_counts)[marker_counts >= min_markers])
drop_clusters <- setdiff(obj_clusters, keep_clusters)

if (length(drop_clusters) > 0) {
  n_drop <- sum(as.character(scobj@meta.data[[cov_markers]]) %in% drop_clusters)
  message("Excluding ", length(drop_clusters), " cluster(s) without >= ", min_markers,
          " markers from CyteTypeR: ", paste(drop_clusters, collapse = ", "),
          " (", n_drop, " cells)")
}
if (length(keep_clusters) == 0) stop("No cluster has >= ", min_markers, " markers; cannot run CyteTypeR.")

keep_cells <- colnames(scobj)[as.character(scobj@meta.data[[cov_markers]]) %in% keep_clusters]
scobj_ct <- subset(scobj, cells = keep_cells)
# re-level so no empty factor levels survive the subset
scobj_ct@meta.data[[cov_markers]] <- droplevels(factor(scobj_ct@meta.data[[cov_markers]]))

df_markers_fix <- df_markers_fix %>% filter(as.character(cluster) %in% keep_clusters)

# run the query preparation
prepped_data <- PrepareCyteTypeR(obj = scobj_ct,
                                 marker_table = df_markers_fix,
                                 n_top_genes = 10,
                                 group_key = cov_markers,
                                 aggregate_metadata = TRUE,
                                 coordinates_key = "umap",
                                 vars_h5_path = out_id_varsh5,
                                 obs_duckdb_path = out_id_duckdb)

# run the annotation
# results <- CyteTypeR(obj = scobj,
#                      auth_token = token,
#                      prepped_data = prepped_data, 
#                      study_context = list_metadata$study_context, 
#                      metadata = list_metadata,
#                      query_filename = out_id_filename)

# always recover output data via GetResults, regardless of whether the direct CyteTypeR() call reports success, since its own return object is missing columns (e.g. granularAnnotation) that the GetResults + transform path always produces. This keeps the annotation table shape consistent across runs. handle the console output to collect the job_id generated
msg_log <- tempfile()
msg_con <- file(msg_log, open = "wt")
sink(msg_con, type = "message")

# submit the query to CyteTypeR; its return value is discarded, it is only used to trigger the job and surface the job_id via the captured log
tryCatch(
  {
    message("Submitting query to CyteTypeR server...")
    CyteTypeR(obj = scobj_ct,
              auth_token = token,
              prepped_data = prepped_data,
              study_context = list_metadata$study_context,
              metadata = list_metadata,
              query_filename = out_id_filename)
  },
  error = function(e) {
    message("CyteTypeR failed: ", conditionMessage(e))
    NULL
  },
  finally = {
    sink(type = "message")
    close(msg_con)
  }
)

# replay captured log so it is printed in the console
logged <- readLines(msg_log)
cat(paste(logged, collapse = "\n"), "\n")

# extract job_id from captured output
all_log <- paste(logged, collapse = " ")
job_id <- regmatches(
  all_log,
  regexpr("[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", all_log)
)
if (length(job_id) == 0) job_id <- NULL
if (is.null(job_id)) {
  stop("CyteTypeR failed before a job_id was obtained. Cannot retrieve results.")
}

# Reconstruct job details in scobj@misc so GetResults uses the same pipeline as CyteTypeR (normalize + correct endpoint), instead of the bare job_id path.
report_url <- paste0("https://prod.cytetype.nygen.io/report/", job_id)
scobj_ct@misc[["cytetype_jobDetails"]] <- list(
  job_id = job_id,
  report_url = report_url,
  api_url = "https://prod.cytetype.nygen.io",
  group_key = cov_markers,
  cluster_labels = prepped_data$clusterLabels
)
message("Retrieving results via GetResults (job_id: ", job_id, ")...")
meta_results <- GetResults(obj = scobj_ct,
                           auth_token = token,
                           api_url = "https://prod.cytetype.nygen.io",
                           results_prefix = "cytetype")

# GetResults returns the normalized result list, not an annotated Seurat object.
# Transform and apply annotations to scobj to mirror what CyteTypeR does internally. cluster_map handles the sequential ID -> original label translation internally.
transformed_results <- CyteTypeR:::.transform_results_seurat(
  meta_results,
  cluster_map = prepped_data$clusterLabels
)

# add the main annotations from the collected results
ann_colname <- paste("cytetype", "annotation",cov_markers, sep = "_")
onto_colname <- paste("cytetype", "cellOntologyTerm", cov_markers, sep = "_")
ontoID_colname <- paste("cytetype", "cellOntologyTermID", cov_markers, sep = "_")
state_colname <- paste("cytetype", "cellState", cov_markers, sep = "_")
granular_colname <- paste("cytetype", "granularAnnotation", cov_markers, sep = "_")

cells_group <- as.character(scobj@meta.data[[cov_markers]])
scobj@meta.data[[ann_colname]] <- factor(setNames(transformed_results$annotation, transformed_results$clusterId)[cells_group])
scobj@meta.data[[onto_colname]] <- factor(setNames(transformed_results$ontologyTerm, transformed_results$clusterId)[cells_group])
scobj@meta.data[[ontoID_colname]] <- factor(setNames(transformed_results$ontologyTermID, transformed_results$clusterId)[cells_group])
scobj@meta.data[[state_colname]] <- factor(setNames(transformed_results$cellState, transformed_results$clusterId)[cells_group])
scobj@meta.data[[granular_colname]] <- factor(setNames(transformed_results$granularAnnotation, transformed_results$clusterId)[cells_group])

# CyteType can return empty fields (e.g. no cellOntologyTerm/ID for a cluster).
# loupeR rejects empty group names downstream, so fill them with "missing".
cytetype_cols <- c(ann_colname, onto_colname, ontoID_colname, state_colname, granular_colname)
scobj@meta.data <- scobj@meta.data %>%
  mutate(across(all_of(cytetype_cols),
                ~ as.character(.x) %>%
                  str_trim() %>%
                  na_if("") %>%
                  replace_na("missing") %>%
                  factor()))

# cells from clusters excluded from the CyteTypeR run (too few markers) get an explicit label instead of the generic "missing"
if (length(drop_clusters) > 0) {
  dropped <- as.character(scobj@meta.data[[cov_markers]]) %in% drop_clusters
  scobj@meta.data <- scobj@meta.data %>%
    mutate(across(all_of(cytetype_cols),
                  ~ factor(replace(as.character(.x), dropped, "Unassigned_no_markers"))))
}

# for some reason in the default run of CyteTypeR there are also some other covariates added to the metadata but are redundant
# cytetype_RNA_snn_res.x -> annoation
# cytetype_ontologyTerm_RNA_snn_res.x -> cellOntologyTerm
# cytetype_ontologyID_RNA_snn_res.x -> cellOntologyTermID

# I have deceded to skip them

scobj@misc$cytetype_results <- transformed_results

results <- scobj

# table(results@meta.data$RNA_snn_res.0.5,results@meta.data$cytetype_cellState_RNA_snn_res.0.5)

# ======================================================================
# == save output ==
# ======================================================================

# save only the per-barcode annotation columns, mirroring AnnotationSingleR's annotation_meta output, so addAnnotations can merge annotations from any combination of methods without needing the full CyteTypeR object.
cytetype_meta <- results@meta.data %>% select(starts_with("cytetype_"))
write.csv(cytetype_meta, out_id_meta)

# save the job_id/report_url separately, since addAnnotations rebuilds the integrated object from scratch and never carries over scobj@misc, so this is the only place downstream steps (e.g. the report) can recover the CyteType job link.
jsonlite::write_json(
  list(job_id = job_id, report_url = report_url, api_url = "https://prod.cytetype.nygen.io"),
  path = out_id_jobinfo,
  auto_unbox = TRUE
)
