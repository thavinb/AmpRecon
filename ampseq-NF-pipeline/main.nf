#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { bcl_to_cram } from './pipeline_workflows/step1.1-bcl-to-cram.nf'
include { get_tag_list_file } from './modules/manifest2tag.nf'

// logging info ----------------------------------------------------------------
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
         --bcl_dir_path            : ${params.bcl_dir_path}

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


// logging info ----------------------------------------------------------------


workflow set_up_params {
    main:
        Channel.fromPath("${params.bcl_dir_path}").set{bcls}
        Channel.value("${params.lane}").set{lane}
        Channel.value("${params.run_id}").set{run_id}
        Channel.value("${params.study_name}").set{study_name}
        Channel.fromPath("${params.manifest_file_path}").set{manifest}

	//manifest2tag
        get_tag_list_file(manifest, "library", "sample", study_name)
        barcodes_file = get_tag_list_file.out.taglist_file

    emit:
        bcls
        lane
        run_id
        study_name
        barcodes_file
}

// Main entry-point workflow

workflow {

    set_up_params()

// Stage 1 - Step 1: BCL to CRAM

    bcl_to_cram(set_up_params.out.bcls, set_up_params.out.lane, set_up_params.out.run_id, set_up_params.out.study_name)//, set_up_params.out.barcodes_file)
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
