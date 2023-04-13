#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows

include { parse_panel_settings } from './modules/parse_panels_settings.nf'
include { SANGER_IRODS_TO_READS } from './workflows/sanger_irods_to_reads.nf'
include { MISEQ_TO_READS } from './workflows/miseq_to_reads.nf'
include { READS_TO_VARIANTS } from './workflows/reads_to_variants.nf'
include { VARIANTS_TO_GRCS } from './workflows/variants_to_grcs.nf'
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
         --execution_mode     : ${params.execution_mode}
         --panels_settings    : ${params.panels_settings}
         --containers_dir     : ${params.containers_dir}
         --results_dir        : ${params.results_dir}
         --containers_dir     : ${params.containers_dir}
         --genotyping_gatk    : ${params.genotyping_gatk}
         --genotyping_bcftools: ${params.genotyping_bcftools}
         --grc_settings_file_path: ${params.grc_settings_file_path}
         --chrom_key_file_path: ${params.chrom_key_file_path}
         --kelch_reference_file_path: ${params.kelch_reference_file_path}
         --codon_key_file_path: ${params.codon_key_file_path}
         --drl_information_file_path: ${params.drl_information_file_path}

         (in-country)
         --run_id             : ${params.run_id}
         --bcl_dir            : ${params.bcl_dir}
         --lane               : ${params.lane}
         --study_name         : ${params.study_name}
         --read_group         : ${params.read_group}
         --library            : ${params.library}

         (irods)
         --irods_manifest     : ${params.irods_manifest}

         (s3)
         --download_from_s3   : ${params.download_from_s3}
         --upload_to_s3       : ${params.upload_to_s3}
         --s3_uuid            : ${params.s3_uuid}
         --s3_bucket_input    : ${params.s3_bucket_input}
         --s3_bucket_output   : ${params.s3_bucket_output}

         (DEBUG)
         --DEBUG_tile_limit   : ${params.DEBUG_tile_limit}
         --DEBUG_takes_n_bams : ${params.DEBUG_takes_n_bams}

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
    (irods)
    nextflow run /path/to/ampseq-pipeline/main.nf -profile sanger_lsf 
      --execution_mode irods
      --irods_manifest ./input/irods_smallset.tsv
      --chrom_key_file_path chromKey.txt
      --grc_settings_file_path grc_settings.json
      --drl_information_file_path DRLinfo.txt
      --codon_key_file_path codonKey.txt
      --kelch_reference_file_path kelchReference.txt
      --containers_dir ./containers_dir/ 

    (incountry)
    nextflow /path/to/ampseq-pipeline/main.nf -profile sanger_lsf
      --execution_mode in-country --run_id 21045
      --bcl_dir /path/to/my_bcl_dir/ --lane 1
      --study_name test --read_group rg_test --library lib
      --chrom_key_file_path chromKey.txt
      --grc_settings_file_path grc_settings.json
      --drl_information_file_path DRLinfo.txt
      --codon_key_file_path codonKey.txt
      --kelch_reference_file_path kelchReference.txt
      --containers_dir ./containers_dir/ 

  Description:
    Ampseq is a bioinformatics analysis pipeline for amplicon sequencing data.
    Currently supporting alignment and SNP variant calling on paired-end Illumina sequencing data.

    *for a complete description of input files and parameters check:
    https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline/

  Options:
    Inputs:
      (required)
      --execution_mode : sets the entry point for the pipeline ("irods" or "in-country")
      
      (incountry required)
      --run_id : id to be used for the batch of data to be processed
      --bcl_dir: path to a miseq directory
      --lane : <str>
      --study_name : <str>
      --read_group : <str>
      --library : <str>

      (irods required)
      --irods_manifest : an tsv containing information of irods data to fetch
      
      (if s3)
      --s3_uuid : <str> a universally unique id which will be used to fetch data from s3, if is not provided, the pipeline will not retrieve miseq runs from s3
      --s3_bucket_input : <str> s3 bucket name to fetch data from
      --upload_to_s3 : <bool> sets if needs to upload output data to an s3 bucket
      --s3_bucket_output : <str> s3 bucket name to upload data to

      (grc_creation)
      --grc_settings_file_path: <str> path to the GRC settings file.
      --chrom_key_file_path: <str> path to the chrom key file
      --kelch_reference_file_path: <str> path to the kelch13 reference sequence file
      --codon_key_file_path: <str> path to the codon key file
      --drl_information_file_path: <str> path to the drug resistance loci information file

    Settings:
      --results_dir : <path>, output directory (Default: $launchDir/output/)
      --panels_settings : <path>, path to panel_settings.csv
      --containers_dir : <path>, path to a dir where the containers are located

      (genotyping)
      --gatk3: <str> path to GATK3 GenomeAnalysisTK.jar file
      --

    Additional options:
      --help (Prints this help message. Default: false)
    
    Profiles:
      sanger_lsf : run the pipeline on farm5 lsf (recommended)
      sanger_default : run the pipeline on farm5 local settings (only for development)
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

  // -- MAIN-EXECUTION ------------------------------------------------------
  // prepare panel resource channels
  ref_and_annt_ch = parse_panel_settings(params.panels_settings)
  reference_ch = ref_and_annt_ch[0] // tuple(reference_file, panel_name, snp_list)
  annotations_ch = ref_and_annt_ch[1] // tuple(panel_name, design_file)

  // Files required for GRC creation
  Channel.fromPath(params.grc_settings_file_path, checkIfExists: true)
  chrom_key_file = Channel.fromPath(params.chrom_key_file_path, checkIfExists: true)
  kelch_reference_file = Channel.fromPath(params.kelch_reference_file_path, checkIfExists: true)
  codon_key_file = Channel.fromPath(params.codon_key_file_path, checkIfExists: true)
  drl_information_file = Channel.fromPath(params.drl_information_file_path, checkIfExists: true)

  if (params.execution_mode == "in-country") {
    // process in country entry point
    MISEQ_TO_READS(reference_ch)
    bam_files_ch = MISEQ_TO_READS.out.bam_files_ch
    file_id_reference_files_ch = MISEQ_TO_READS.out.file_id_reference_files_ch
    file_id_to_sample_id_ch = MISEQ_TO_READS.out.file_id_to_sample_id_ch
  }

  if (params.execution_mode == "irods") {
    // process IRODS entry point
    SANGER_IRODS_TO_READS(params.irods_manifest, reference_ch)
    // setup channels for downstream processing
    bam_files_ch = SANGER_IRODS_TO_READS.out.bam_files_ch // tuple (file_id, bam_file, run_id)
    file_id_reference_files_ch = SANGER_IRODS_TO_READS.out.file_id_reference_files_ch
    file_id_to_sample_id_ch = SANGER_IRODS_TO_READS.out.file_id_to_sample_id_ch
  }

  // Reads to variants
  READS_TO_VARIANTS(bam_files_ch, file_id_reference_files_ch, annotations_ch, file_id_to_sample_id_ch)
  lanelet_manifest_file = READS_TO_VARIANTS.out.lanelet_manifest

  // Variants to GRCs
  VARIANTS_TO_GRCS(lanelet_manifest_file, chrom_key_file, kelch_reference_file, codon_key_file, drl_information_file)

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

def validate_general_params(){
  /*
  count errors on parameters which must be provided regardless of the workflow which will be executed
  
  returns
  -------
  <int> number of errors found
  */

  def err = 0
  def valid_execution_modes = ["in-country", "irods"]

  // check if execution mode is valid
  if (!valid_execution_modes.contains(params.execution_mode)){
    log.error("The execution mode provided (${params.execution_mode}) is not valid. valid modes = ${valid_execution_modes}")
    err += 1
  }
  // check if output dir exists, if not create the default
  if (params.results_dir){
    results_path = file(params.results_dir)
    if (!results_path.exists()){
      log.warn("${results_path} does not exists, the dir will be created")
      results_path.mkdir()
    }
  }
  
  // if S3 is requested, check if all s3 required parameters were provided
  // check S3 input bucket
  if (!(params.s3_bucket_input==null)){
    if (params.s3_uuid == null){
      log.error("A s3 uuid parameter must be provided if a s3 bucket input is provided'.")
      errors += 1
    } 
  }
  // check S3 output bucket
  if (params.upload_to_s3){
    if (params.s3_bucket_output == null){
      log.error("A s3_bucket_output parameter must be provided if upload_to_s3 is set to '${params.upload_to_s3}'.")
      errors += 1
    }
    if (params.s3_uuid==null){
      log.error("A s3_uuid must be provided if upload_to_s3 is required")
      errors += 1
    }
  }
  return err
}
