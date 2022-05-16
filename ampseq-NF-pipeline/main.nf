#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { reads_retrieval_workflow  } from './pipeline_workflows/reads_retrieval_workflow.nf'
include { filtering_workflow } from './pipeline_workflows/filtering_workflow.nf'
include { variant_calling_workflow } from './pipeline_workflows/variant_calling_workflow.nf'
include { genotyping_workflow } from './pipeline_workflows/genotyping_workflow.nf'

// logginh info ----------------------------------------------------------------
// This part of the code is based on the FASTQC PIPELINE (https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf)

/*
* ANSI escape codes to color output messages, get date to use in results folder name
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

log.info """
        ===========================================
         AMPSEQ_0.0 (dev : prototype)
         Used parameters:
        -------------------------------------------
         --bcl_path            : ${params.bcl_path}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Run container:          ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
         Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         ------------------------------------------
         """
         .stripIndent()


// logginh info ----------------------------------------------------------------


workflow set_up_params {
    main:
        Channel.fromPath("${params.bcl_path}/*.bcl").set {bcls}
        //Channel.fromPath(params.irods_paths).set {irods_paths}
    emit:
        bcls
        //irods_paths
}

// Main entry-point workflow

workflow {

    set_up_params()

// Step 1 - Reads retrieval

    reads_retrieval_workflow(set_up_params.out)

// Step 2 - Added filtering

    filtering_workflow(reads_retrieval_workflow.out)

// Step 3 - Variant calling

    variant_calling_workflow(filtering_workflow.out)

// Step 4 - Drug resistance genotyping

    genotyping_workflow(variant_calling_workflow.out)
}

// -------------- Check if everything went okay -------------------------------
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            ${ANSI_GREEN}Finished in ${workflow.duration}
            See the report here ==> ${ANSI_RESET}/SOMEDIR/XXX_report.html
            """
            .stripIndent()
    } else {
        log.info """
            ===========================================
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}
