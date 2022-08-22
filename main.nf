#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows

include { IN_COUNTRY } from './workflows/in_country.nf'

// - process to extract and validate information expected based on input params
include { validate_parameters } from './workflows/pipeline-subworkflows/inputHandling.nf'


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
         --input_params_csv           : ${params.input_params_csv}
         ------------------------------------------
         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         ------------------------------------------
         """
         .stripIndent()


def printHelp() {
  log.info """
  Usage:
    nextflow run main.nf --input_params_csv [path/to/my/input.csv]

  Description:
    (temporary - honest - description)
    Ampseq amazing brand new pipeline built from the ground up to be awesome.

    A input csv containing my run_id, bcl_dir, study_name, read_group, library,
    and reference fasta path is necessary to run the ampseq pipeline from step 0

    *for a complete description of input files check [LINK AMPSEQ]

  Options:
    Inputs:
      --input_params_csv (A csv file path)
      --irods_manifest (tab-delimited file containing rows of WG_lane and id_run data for CRAM files on iRODS)

    Additional options:
      --help (Prints this help message. Default: false)
      --results_dir (Results directory. Default: $launchDir/output/)
   """.stripIndent()
}

// Main entry-point workflow
workflow {
  // --- Print help if requested -------------------------------------------
  // Show help message
  if (params.help) {
      printHelp()
      exit 0
  }

  // -- MAIN-EXECUTION -------------------------------------------------------------

  if (params.mode == "in-country") {
        IN_COUNTRY()
  }

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
