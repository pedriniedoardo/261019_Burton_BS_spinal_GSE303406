rule sampleMetadata:
  '''
  Generate a simple metadata file (from samples.csv). This file is meant to
  be hand-edited afterwards with the real sample annotations, so it is
  intentionally NOT wired to LUT_QC: which samples actually end up in the
  report is decided downstream, in reportImages, by filtering against
  LUT_QC's include_sample column (see 05_report_images.R). Keeping this
  rule free of that dependency means editing LUT_QC never triggers Snakemake
  to regenerate (and wipe) your hand-typed annotations.
  '''
  params:
    config["report"]["sample_annotations"]
  output:
    config["samples_meta_csv"]
  run:
    # read samples config
    df_sample = pd.read_csv(config["samples_csv"])
    # base metadata sempre derivabile
    df_base_meta = (
        df_sample[['sample_name', 'sample_run']]
        .drop_duplicates()
        .groupby('sample_name')['sample_run']
        .agg(';'.join)
        .reset_index()
    )
    # add constant annotation columns
    for col in params[0]:
        df_base_meta[col] = "all"
    # write to disk
    df_base_meta.to_csv(output[0], index=False)
    
rule reportImages:
  '''
  This rule is intended for producing images on the embedded Seurat's object
  '''
  input:
    pre_qc_rds     = expand(rules.runSamplePreprocessing.output.rds,
          sample_name=SAMPLES_merge.keys()),
    merged_rds     = rules.integrateSamples.output.merged_rds,
    integrated_rds = rules.addAnnotations.output.rds,
    samples_meta   = rules.sampleMetadata.output,
    LUT_QC         = rules.runAggregateQC.output.LUT_QC,
    markers_tsv    = rules.integrateSamples.output.markers
  output:
    tab_qc = config["out_location"] + config["out_location_seurat"] +
        'table/stats_preQC_postQC_' + CELLBENDER_TAG + '.csv',
    done = config["out_location"] + config["out_location_seurat"] + "report_" + CELLBENDER_TAG + ".done"
    # # ----------- table of the quality controls ----------- #
    # tab_qc = config["out_location"] + 
    #     'Seurat/table/stats_preQC_postQC.csv',
    # # ----------- violin plots ---------------------------- #
    # vln_qc = expand(
    #     config["out_location"] + config["out_location_seurat"] + "plot/vln_" + CELLBENDER_TAG + ".{extension}",
    #     extension = config["report"]["figure_extension"]
    # ),
    # # ----------- UAMP with quality chekcings ------------- #
    # dim_qc = expand(
    #     config["out_location"] + config["out_location_seurat"] + "plot/dim_{feature}_" + CELLBENDER_TAG + ".{extension}",
    #     feature   = config["report"]["qc_features"],
    #     extension = config["report"]["figure_extension"]
    # ),
    # # ----------- UAMP with samples colored --------------- #
    # dim_ident = expand(
    #   config["out_location"] + config["out_location_seurat"] + "plot/dim_{object_type}_ident_" + CELLBENDER_TAG + ".{extension}",
    #   object_type = ["merged", "integrated"],
    #   extension   = config["report"]["figure_extension"]
    #   ),
    # # ----------- UAMP with cluster colored --------------- #
    # dim_clusters = expand(
    #   config["out_location"] + config["out_location_seurat"] + "plot/dim_{object_type}_cluster_" + CELLBENDER_TAG + ".{extension}",
    #   object_type = ["merged", "integrated"],
    #   resolution  = config["Seurat"]["FindClusters"]["resolution"],
    #   extension   = config["report"]["figure_extension"]
    #   ),
    # # ----------- UAMP with cluster colored, split by ann - #
    # # dim_clusters_byann = expand(
    # #   config["out_location"] + 
    # #       "Seurat/plot/dim_{object_type}_cluster_RES{resolution}_by{ann}_" + 
    # #       CELLBENDER_TAG + ".{extension}",
    # #   object_type = ["merged", "integrated"],
    # #   resolution  = config["Seurat"]["FindClusters"]["resolution"],
    # #   ann         = config["report"]["sample_annotations"],
    # #   extension   = config["report"]["figure_extension"]
    # #   ),
    #   # ----------- UAMP with cell.cycle colored ---------- #
    # dim_cellcycle = expand(
    #   config["out_location"] + config["out_location_seurat"] + "plot/dim_{object_type}_CellCycle_" + 
    #            CELLBENDER_TAG + ".{extension}",
    #   object_type = ["merged", "integrated"],
    #   extension   = config["report"]["figure_extension"]
    #   )
  params:
    out_location     = config["out_location"] + config["out_location_seurat"],
    qc_features      = config["report"]["qc_features"],
    figure_extension = config["report"]["figure_extension"],
    # ann_cols         = config["report"]["sample_annotations"],
    resolutions      = config["Seurat"]["FindClusters"]["resolution"],
    cellbender_tag   = CELLBENDER_TAG
  conda: config["env_seurat"]
  log:
      'logs/report/reportImages.log'
  benchmark:
      'benchmarks/report/reportImages.txt'
  script:
      "../scripts/05_report_images.R"
