rule renderReport:
    """
    Render the scRNAseq HTML report from the R Markdown template.
    Depends on all figure and table outputs produced by reportImages.
    """
    input:
        seurat_rds = rules.addAnnotations.output.rds,
        rmd    = "report/report.Rmd",
        cytetype_job_id_file = rules.AnnotationCyteTypeR.output.job_id_file,
        samples_meta = rules.sampleMetadata.output,
        tab_qc = rules.reportImages.output.tab_qc,
        done   = rules.reportImages.output.done,
    output:
        html = "report/report_" + CELLBENDER_TAG + ".html"
    conda:
        config["env_seurat"]
    params:
        author = config["report"]["author"],   # "Aurora Maurizio - Center for Omics Sciences (COSR)"
        genome = config["ref"]["build"],
        nPCs       = config["Seurat"]["FindNeighbors"]["dims"],
        resolutions = lambda wildcards: "c(" + ",".join(str(r) for r in config["Seurat"]["FindClusters"]["resolution"]) + ")",
        report_resolution = config["report"]["report_resolution"],
        cellbender_tag = CELLBENDER_TAG,
        seurat_rds = rules.addAnnotations.output.rds,
        out_location_seurat = config["out_location"] + config["out_location_seurat"]
    retries: 2  # re-run on failure (e.g. OOM); mem_mb in the profile scales with `attempt`
    log:
        "logs/renderReport/renderReport.log"
    benchmark:
        "benchmarks/renderReport/renderReport.txt"
    shell:
        """
        cd report && \
        PANDOC_STACK_SIZE=2048m \
        LC_ALL=C.UTF-8 \
        Rscript -e "rmarkdown::render('report.Rmd', \
          params = list(
            author     = '{params.author}',
            genome     = '{params.genome}',
            nPCs       = {params.nPCs},
            resolutions = {params.resolutions},
            report_resolution = {params.report_resolution},
            cellbender_tag = '{params.cellbender_tag}',
            out_location_seurat = '../{params.out_location_seurat}',
            cytetype_job_id_file = '../{input.cytetype_job_id_file}',
            samples_meta = '../{input.samples_meta}',
            seurat_rds = '../{params.seurat_rds}' \
          ), \
          output_file = 'report_{params.cellbender_tag}.html')" \
        > "../{log}" 2>&1
        """
