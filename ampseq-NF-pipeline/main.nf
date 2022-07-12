#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows
include { prepare_reference } from './pipeline_workflows/step0.1b-prepare-reference.nf'
include { bcl_to_cram } from './pipeline_workflows/step1.1-bcl-to-cram.nf'
include { cram_to_bam } from './pipeline_workflows/step1.2-cram-to-bam.nf'
include { reset_bam_alignment } from './pipeline_workflows/step1.3-bam-to-vcf'
include { pull_from_iRODS } from './pipeline_workflows/step1.2b-pull-from-iRODS.nf'

// - process to extract and validate information expected based on input params
include { get_taglist_file } from './modules/manifest2tag.nf'
include { validate_parameters; load_manifest_ch; load_steps_to_run } from './modules/inputHandling.nf'
include { make_samplesheet_manifest } from './modules/make_samplesheet_manifest.nf'
include { validate_samplesheet_manifest } from './modules/samplesheet_manifest_validation.nf'
include { samplesheet_validation } from './modules/samplesheet_validation.nf'

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
         --manifest        (required) : ${params.manifest}
         --reference_fasta (required) : ${params.reference_fasta}
         --start_from                 : ${params.start_from}
         --results_dir                : ${params.results_dir}
         --irods_manifest             : ${params.irods_manifest}

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
    nextflow run main.nf --manifest [path/to/my/manifest.csv]
    [TO DO] nextflow run main.ng --from_step 1.2 --manifest_step1_2 [path/to/my/manifest_for_step1_2.csv]

  Description:
    (temporary - honest - description)
    Ampseq amazing brand new pipeline built from the ground up to be awesome.

    A manifest containing my run_id, bcl_dir, study_name, read_group, library,
    and reference fasta path is necessary to run the ampseq pipeline from step 0

    *for a complete description of input files check [LINK AMPSEQ]

  Options:
    Inputs:
      --manifest (A manifest csv file)
      --manifest_step1_2 (manifest file to be submitted to step 1.2, previous steps are ignored)
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

  // -- Pre Processing ------------------------------------------------------
  // validate input
  validate_parameters()
  // get steps to run
  steps_to_run_tags = load_steps_to_run()
  tag_provided = params.start_from.toString()
  println(steps_to_run_tags)
  if (steps_to_run_tags.contains("0")) {
    // process manifest input
    manifest_ch = load_manifest_ch()
    // validate MiSeq run files and directory structure
    samplesheet_validation(manifest_ch)

    // process samplesheets manifest (necessary to get barcodes) and validate it
    make_samplesheet_manifest(manifest_ch)//run_id, manifest_ch.bcl_dir)
    validate_samplesheet_manifest(make_samplesheet_manifest.out)
    // get taglist
    get_taglist_file_In_ch = manifest_ch.join(validate_samplesheet_manifest.out)
    get_taglist_file(get_taglist_file_In_ch)

    // set step1 input channel
    step1_Input_ch = manifest_ch.join(get_taglist_file.out)

    // Stage 1 - Step 1: BCL to CRAM
    bcl_to_cram(step1_Input_ch)
    manifest_step1_1_Out_ch = bcl_to_cram.out.multiMap { it -> run_id: it[0]
                                                                  mnf: it[1]}
  }
  if (steps_to_run_tags.contains("1.2a")) {
    // handle reference fasta
    prepare_reference(params.reference_fasta)
    reference_idx_fls = prepare_reference.out
  }

  if (steps_to_run_tags.contains("1.2a")) {

    // if start from this step, use the provided in_csv, if not, use previous
    // step output
    if (tag_provided=="1.2a"){
      step1_2_In_mnf = params.step1_2_in_csv
      csv_ch = Channel.fromPath(params.step1_2_in_csv)
    }
    else {
      csv_ch = manifest_step1_1_Out_ch.mnf
    }

    // Stage 1 - Step 2: CRAM to BAM
    cram_to_bam(csv_ch,
                reference_idx_fls.bwa_index_fls,
                reference_idx_fls.fasta_index_fl,
                reference_idx_fls.dict_fl)//,
                //irods_ch)
    cram_to_bam.out
    step1_2_Out_ch = cram_to_bam.out.multiMap { it -> sample_tag: it[0]
                                                      bam_file: it[1]
                                                      run_id:it[2]}
  }
   // --- ADD irods pulling here 1.2b ---------------------------------------
   // it should generate an 1.2b_out_csv equivalent to the 1.2a_out_csv
   // -----------------------------------------------------------------------
   //

   if (tag_provided=="1.2b"){

     irods_ch =  Channel.fromPath(params.irods_manifest, checkIfExists: true)
                      | splitCsv(header: true, sep: '\t')
                      | map { row -> tuple(row.id_run, row.WG_lane) }

     // run step1.2b - pull from iRODS
     pull_from_iRODS(irods_ch)
     // prepare channel for step 1.3
     step1_2_Out_ch = pull_from_iRODS.out.multiMap { it -> sample_tag: it[0]
                                                           bam_file: it[1]
                                                           run_id:it[2]
                                                   }
   }


  if (steps_to_run_tags.contains("1.3")){
    // if start from this step, use the provided in_csv, if not, use step 1.2x
    // step output
    if (tag_provided=="1.3"){
      step1_3_In_mnf = params.step1_3_in_csv
      step1_3_In_ch = step1_3_In_mnf.splitCsv(header : true)
                            .multiMap {
                              row  -> run_id:row.run_id
                                      bam_file:row.bam_fl
                                      sample_tag:row.sample_tag
                }
    } else {
        step1_3_In_ch = step1_2_Out_ch
    }
    // run step1.3 - BAM to VCF
    reset_bam_alignment(step1_3_In_ch.sample_tag,
                        step1_3_In_ch.bam_file,
                        step1_3_In_ch.run_id)
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
