process Banksy {
    publishDir "${params.outdir}/",
        mode: 'copy'
    tag "markers"
    label "highMemMT1"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == "apptainer" ?
        'docker://yuanli202004/seurat5.4_doubletfinder:v1.0.3':
        'docker://yuanli202004/seurat5.4_doubletfinder:v1.0.3' }"

    input:
        path(INFILE)

    output:
        path("celltypeNei"), emit: outDir
        
    script:
    """
    mkdir -p celltypeNei
    Rscript ${moduleDir}/templates/clusterNeighborhoodAnalysis.R \
        --OUTDIR "celltypeNei/" \
        --INFILE "${INFILE}"

    """
}

