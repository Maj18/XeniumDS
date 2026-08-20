#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    XeniumDS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : 
    Contact: 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

println """\
         XeniumDS   P I P E L I N E
         ===================================
         GitHub: 
         ___________________________________
         OUTPUT DIR     : ${params.outdir}
        ___________________________________
         """
         .stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { getDataCSV } from './modules/getDataCSV'
include { clusterNei } from './modules/clusterNeighborhoodAnalysis'

workflow ENTRY_clusterNei {
    def INFILE = Channel.value(params.INFILE)
    clusterNei(INFILE)
}

workflow {
    // **************** 1. Save layer normalized data to scv files ****************
    ch_INFILE = Channel.value(params.INFILE)
    getDataCSV(ch_INFILE)
    ch_layer_csvs_dir = getDataCSV.out.layer_csvs_dir

}

workflow.onComplete {
    println( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

