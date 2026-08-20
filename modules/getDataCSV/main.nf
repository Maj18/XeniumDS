process getDataCSV {
    publishDir "${params.outdir}/data/",
        mode: 'copy'
    tag "data"
    label "highMemMT1"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'docker://yuanli202004/seurat5.4_doubletfinder:latest':
        'docker://yuanli202004/seurat5.4_doubletfinder:latest' }"

    input:
        path(INFILE)

    output:
        path("layer_csvs/"), emit: layer_csvs_dir

    script:
    """
    Rscript ${moduleDir}/templates/QC.R \
        --INFILE "${INFILE}"

    """
}

