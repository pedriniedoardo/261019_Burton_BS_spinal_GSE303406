rule runSamplePreprocessing:
    '''
    This rule preprcess the individual samples from the cellranger output.
    '''
    input:
        # result_folder = rules.runCellrangerMultiRun.output.folder
        target = get_preprocessing_input
    output:
        rds = config["out_location"] + config["out_location_seurat"] + "object/{sample_name}_" + CELLBENDER_TAG + "_obj_preQC.rds",
        meta = config["out_location"] + config["out_location_seurat"] + "table/{sample_name}_" + CELLBENDER_TAG + "_meta_preQC.tsv",
    conda: config["env_seurat"]
    log:
        'logs/Seurat/{sample_name}/runSamplePreprocessing.log'
    benchmark:
        'benchmarks/Seurat/{sample_name}/runSamplePreprocessing.txt'
    params:
        id_org = config["ref"]["organism_id"],
        # used to make update the implementation of reading the h5 file.
        # Pass the boolean flag to R
        use_cellbender = config.get("run_cellbender", False)
    script:
        "../scripts/01_sample_preprocessing_snakemake.R"

rule runAggregateQC:
    '''
    This rule aggreagates the metadata form the individual samples and generates the plot to define the QC thresholds.
    '''
    input:
        # this is needed to trigger it after the generation of the outputs on all the outputs
        meta_tables = expand(rules.runSamplePreprocessing.output.meta,sample_name=SAMPLES_merge.keys())
    output:
        plot_mito = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_mito_V5_" + CELLBENDER_TAG + ".pdf",
        plot_ribo = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_ribo_V5_" + CELLBENDER_TAG + ".pdf",
        plot_hemo = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_hemo_V5_" + CELLBENDER_TAG + ".pdf",
        plot_count = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_count_V5_" + CELLBENDER_TAG + ".pdf",
        plot_feature = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_features_V5_" + CELLBENDER_TAG + ".pdf",
        LUT_QC = config["out_location"] + config["out_location_seurat"] + "table/LUT_QC_" + CELLBENDER_TAG + ".csv",
        meta_total = config["out_location"] + config["out_location_seurat"] + "table/meta_total_beforeQC_V5_" + CELLBENDER_TAG + ".tsv"
    conda: config["env_seurat"]
    log:
        'logs/Seurat/runAggregateQC.log'
    benchmark:
        'benchmarks/Seurat/runAggregateQC.txt'
    params:
        # id_org = config["ref"]["organism_id"]
        # pass all the sample names in the script
        sample_names = list(SAMPLES_merge.keys())
    script:
        "../scripts/01_aggregate_QC_snakemake.R"

rule plotQCThresholds:
    '''
    This rule regenerates the QC diagnostic plots from runAggregateQC, overlaying the per-sample thresholds defined in the (user-edited) LUT_QC table so they can be reviewed before runApplyQC.
    '''
    input:
        meta_total = rules.runAggregateQC.output.meta_total,
        LUT_QC = rules.runAggregateQC.output.LUT_QC
    output:
        plot_mito = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_mito_thr_V5_filters_" + CELLBENDER_TAG + ".pdf",
        plot_ribo = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_ribo_thr_V5_filters_" + CELLBENDER_TAG + ".pdf",
        plot_hemo = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_hemo_thr_V5_filters_" + CELLBENDER_TAG + ".pdf",
        plot_count = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_count_thr_V5_filters_" + CELLBENDER_TAG + ".pdf",
        plot_feature = config["out_location"] + config["out_location_seurat"] + "plot/fixed_histo_features_thr_V5_filters_" + CELLBENDER_TAG + ".pdf"
    conda: config["env_seurat"]
    log:
        'logs/Seurat/plotQCThresholds.log'
    benchmark:
        'benchmarks/Seurat/plotQCThresholds.txt'
    script:
        "../scripts/01_plot_QC_thresholds_snakemake.R"

# checkpoint runApplyQC:
#     '''
#     This rule preprcess the individual samples from the cellranger output.
#     '''
#     input:
#         rds = rules.runSamplePreprocessing.output.rds,
#         LUT_QC = rules.runAggregateQC.output.LUT_QC
#     output:
#         rds = config["out_location"] + config["out_location_seurat"] + "object/{sample_name}_obj_postQC.rds",
#         meta = config["out_location"] + config["out_location_seurat"] + "table/{sample_name}_meta_postQC.tsv",
#     conda: config["env_seurat"]
#     log:
#         'logs/Seurat/{sample_name}/runApplyQC.log'
#     benchmark:
#         'benchmarks/Seurat/{sample_name}/runApplyQC.txt'
#     params:
#         # id_org = config["ref"]["organism_id"]
#     script:
#         "../scripts/02_apply_QC_snakemake.R"

rule runApplyQC:
    '''
    This rule preprcess the individual samples from the cellranger output.
    '''
    input:
        rds = rules.runSamplePreprocessing.output.rds,
        LUT_QC = rules.runAggregateQC.output.LUT_QC
    output:
        rds = config["out_location"] + config["out_location_seurat"] + "object/{sample_name}_" + CELLBENDER_TAG + "_obj_postQC.rds",
        meta = config["out_location"] + config["out_location_seurat"] + "table/{sample_name}_" + CELLBENDER_TAG + "_meta_postQC.tsv",
    conda: config["env_seurat"]
    log:
        'logs/Seurat/{sample_name}/runApplyQC.log'
    benchmark:
        'benchmarks/Seurat/{sample_name}/runApplyQC.txt'
    params:
        # id_org = config["ref"]["organism_id"]
    script:
        "../scripts/02_apply_QC_snakemake.R"


# ====================================================================
# == This function computes the number of cells in each experiment,
# == in order to assign resources for the integration step
# ====================================================================

# def compute_memory_usage(wildcards):
    
#     mem_per_cell = config["memory_per_cell"]
#     max_mem      = config["max_mem_mb"]
    
#     total_cells = 0
#     for sample in SAMPLES_merge.keys():
#         ck = checkpoints.runApplyQC.get(sample_name=sample)
#         meta_file = ck.output["meta"]
#         with open(meta_file) as fp:
#             total_cells += sum(1 for _ in fp)

#     requested_mem = total_cells * mem_per_cell
#     return min(requested_mem, max_mem)

# rule integrateSamples:
#     '''
#     This rule integrates the different samples via Harmony
#     '''
#     input:
#         rds = lambda wc: [checkpoints.runApplyQC.get(sample_name=s).output.rds
#                   for s in SAMPLES_merge.keys()]
#     output:
#         rds = config["out_location"] + config["out_location_seurat"] + "object/integrated_obj.rds",
#         meta = config["out_location"] + config["out_location_seurat"] + "table/integrated_obj_meta.tsv",
#         markers = config["out_location"] + config["out_location_seurat"] + "table/integrated_obj_markers.tsv",
#     conda: config["env_seurat"]
#     resources:
#         mem_mb = compute_memory_usage
#     log:
#         'logs/Seurat/integrated/runApplyQC.log'
#     benchmark:
#         'benchmarks/Seurat/integrated/integrateSamples.txt'
#     params:
#         # id_org = config["ref"]["organism_id"]
#     script:
#         "../scripts/03_integrate.R"

rule integrateSamples:
    '''
    This rule integrates the different samples via Harmony
    '''
    input:
        rds = expand(rules.runApplyQC.output.rds,
          sample_name=SAMPLES_merge.keys()),
        meta_tables = expand(rules.runApplyQC.output.meta,
          sample_name=SAMPLES_merge.keys()),
        LUT_QC = rules.runAggregateQC.output.LUT_QC
        # samples_meta = "metadata/metadata.csv"
    output:
        merged_rds     = config["out_location"] + config["out_location_seurat"] + "object/merged_obj_" + CELLBENDER_TAG + ".rds",
        integrated_rds = config["out_location"] + config["out_location_seurat"] + "object/integrated_obj_" + CELLBENDER_TAG + ".rds",
        meta           = config["out_location"] + config["out_location_seurat"] + "table/integrated_meta_" + CELLBENDER_TAG + ".tsv",
        markers        = config["out_location"] + config["out_location_seurat"] + "table/integrated_markers_" + CELLBENDER_TAG + ".tsv",
    conda: config["env_seurat"]
    resources:
        mem_mb = lambda wildcards, input: min(
            config["max_mem_mb"],
            sum(sum(1 for line in open(f)) for f in input.meta_tables) * config["memory_per_cell"]
        ) if all(os.path.exists(f) for f in input.meta_tables) else config["max_mem_mb"]
    log:
        'logs/Seurat/integrated/integrateSamples.log'
    benchmark:
        'benchmarks/Seurat/integrated/integrateSamples.txt'
    params:
        seurat_params = config["Seurat"],
        organism_id   = config["ref"]["organism_id"],
        df_cc         = config["df_cc"]
    script:
        "../scripts/03_integrate.R"
        # "../scripts/03_integrate_cellcycle.R"
        # "../scripts/03_integrate_cellcycle_filter.R"

rule AnnotationSingleR:
    '''
    This rule run the automatic annotation using SingleR
    '''
    input:
        rds = rules.integrateSamples.output.integrated_rds,
        meta_tables = expand(rules.runApplyQC.output.meta,
          sample_name=SAMPLES_merge.keys()) # for dynamical allocation of resources
    output:
        annotation_meta = config["out_location"] + config["out_location_seurat"] + "table/integrated_singleR_meta_" + CELLBENDER_TAG + ".csv",
        singleR_output_main  = config["out_location"] + config["out_location_seurat"] + "object/integrated_singleR_main_prediction_" + CELLBENDER_TAG + ".rds",
        singleR_heatmap_main = config["out_location"] + config["out_location_seurat"] + "plot/integrated_singleR_main_plotScoreHeatmap_" + CELLBENDER_TAG + ".png",
        singleR_output_fine  = config["out_location"] + config["out_location_seurat"] + "object/integrated_singleR_fine_prediction_" + CELLBENDER_TAG + ".rds",
        singleR_heatmap_fine = config["out_location"] + config["out_location_seurat"] + "plot/integrated_singleR_fine_plotScoreHeatmap_" + CELLBENDER_TAG + ".png"
    conda: config["env_singleR"]
    resources:
        mem_mb = lambda wildcards, input: min(
            config["max_mem_mb"],
            sum(sum(1 for line in open(f)) for f in input.meta_tables) * config["memory_per_cell"]
        ) if all(os.path.exists(f) for f in input.meta_tables) else config["max_mem_mb"]
    log:
        'logs/Seurat/integrated/AnnotationSingleR.log'
    benchmark:
        'benchmarks/Seurat/integrated/AnnotationSingleR.txt'
    params:
        reference_dataset = config["SingleR_reference"]
    script:
        "../scripts/04_A_annotation_singleR.R"

# token for the annotation step is defined as environment variable
envvars:
        "CYTETYPER_TOKEN"

rule AnnotationCyteTypeR:
    '''
    This rule run the automatic annotation using CyteTypeR
    '''
    input:
        rds = rules.integrateSamples.output.integrated_rds,
        markers = rules.integrateSamples.output.markers,
        metadata = config["meta_annotation"],
        meta_tables = expand(rules.runApplyQC.output.meta,
        sample_name=SAMPLES_merge.keys())
    output:
        annotation_meta = config["out_location"] + config["out_location_seurat"] + "table/integrated_cytetypeR_meta_" + CELLBENDER_TAG + ".csv",
        vars_h5_path = config["out_location"] + config["out_location_seurat"] + "object/vars_" + CELLBENDER_TAG + ".h5",
        obs_duckdb_path = config["out_location"] + config["out_location_seurat"] + "object/obs_" + CELLBENDER_TAG + ".duckdb",
        query_filename = config["out_location"] + config["out_location_seurat"] + "object/query_" + CELLBENDER_TAG + ".json",
        job_id_file = config["out_location"] + config["out_location_seurat"] + "table/cytetype_job_id_" + CELLBENDER_TAG + ".json"
    conda: config["env_annotationCyteTypeR"]
    resources:
        mem_mb = lambda wildcards, input: min(
            config["max_mem_mb"],
            sum(sum(1 for line in open(f)) for f in input.meta_tables) * config["memory_per_cell"]
        ) if all(os.path.exists(f) for f in input.meta_tables) else config["max_mem_mb"]
    log:
        'logs/Seurat/integrated/AnnotationCyteTypeR.log'
    benchmark:
        'benchmarks/Seurat/integrated/AnnotationCyteTypeR.txt'
    params:
        cellbender_tag   = CELLBENDER_TAG,
        # id_org = config["ref"]["organism_id"]
        # token = config["token_CyteTypeR"],
        cov_markers = config["cov_markers_CyteTypeR"]
    script:
        "../scripts/04_annotation_CyteTypeR.R"

# ====================================================================
# == select which automatic annotation method(s) feed the final object
# ====================================================================

ANNOTATION_META_LOOKUP = {
    "singleR":   rules.AnnotationSingleR.output.annotation_meta,
    "cytetypeR": rules.AnnotationCyteTypeR.output.annotation_meta,
}

def get_annotation_metas(wildcards):
    methods = config["annotation_methods"]
    unknown = set(methods) - set(ANNOTATION_META_LOOKUP)
    if unknown:
        raise ValueError(
            f"Unknown annotation_methods {unknown}, choose from {list(ANNOTATION_META_LOOKUP)}"
        )
    return [ANNOTATION_META_LOOKUP[m] for m in methods]

rule addAnnotations:
  '''
  This rule adds the cell type annotations to the integrated object
  '''
  input:
    rds = rules.integrateSamples.output.integrated_rds,
    annotation_metas = get_annotation_metas,
    meta_tables = expand(rules.runApplyQC.output.meta,
      sample_name=SAMPLES_merge.keys()) # for dynamical allocation of resources
  output:
    rds = config["out_location"] + config["out_location_seurat"] + 
      "object/integrated_obj_" + CELLBENDER_TAG + "_annotated.rds",
  conda: config["env_seurat"]
  resources:
      mem_mb = lambda wildcards, input: min(
          config["max_mem_mb"],
          sum(sum(1 for line in open(f)) for f in input.meta_tables) * config["memory_per_cell"]
      ) if all(os.path.exists(f) for f in input.meta_tables) else config["max_mem_mb"]
  log: 
    'logs/Seurat/integrated/addAnnotations.log'
  benchmark: 
    'benchmarks/Seurat/integrated/addAnnotations.txt'
  script: 
    "../scripts/04_add_annotations.R"
    
