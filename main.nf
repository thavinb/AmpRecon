#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows

include { PARSE_PANEL_SETTINGS } from './workflows/parse_panels_settings.nf'
include { IRODS } from './workflows/irods.nf'
include { IN_COUNTRY } from './workflows/in_country.nf'
include { COMMON } from './workflows/common.nf'
include { validate_parameters } from './workflows/input_handling.nf'
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
         --genotyping_gatk    : ${params.genotyping_gatk}
         --genotyping_bcftools: ${params.genotyping_bcftools}
         --panels_settings    : ${params.panels_settings}
         --containers_dir     : ${params.containers_dir}
         --results_dir        : ${params.results_dir}

         (in-country)
         --run_id             : ${params.run_id}
         --bcl_dir            : ${params.bcl_dir}
         --lane               : ${params.lane}
         --study_name         : ${params.study_name}
         --read_group         : ${params.read_group}
         --library            : ${params.library}

         (irods)
         --irods_manifest     : ${params.irods_manifest}
         --aligned_bams_mnf   : ${params.aligned_bams_mnf}
         
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

    (incountry)
    nextflow /path/to/ampseq-pipeline/main.nf -profile sanger_lsf
                --execution_mode in-country --run_id 21045
                --bcl_dir /path/to/my_bcl_dir/ --lane 1
                --study_name test --read_group rg_test --library lib

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

    Settings:
      --results_dir : <path>, output directory (Default: $launchDir/output/)
      --panels_settings : <path>, path to panel_settings.csv
      --containers_dir : <path>, path to a dir where the containers are located

      (genotyping)
      --gatk3: <str> path to GATK3 GenomeAnalysisTK.jar file.
      --combined_vcf_file1 : <path> known SNPs database file. Used to prevent BaseRecalibrator from using regions surrounding polymorphisms.
      --combined_vcf_file2 : <path> known SNPs database file. Used to prevent BaseRecalibrator from using regions surrounding polymorphisms.
      --combined_vcf_file3 : <path> known SNPs database file. Used to prevent BaseRecalibrator from using regions surrounding polymorphisms.
      --conserved_bed_file : <path> file containing genomic intervals the GATK BaseRecalibrator command operates over in the bqsr.nf process.
      --gatk_base_recalibrator_options : <str> input settings containing the supplied known sites files paths and intervals file path for the BaseRecalibrator command in the bqsr.nf process.
      --alleles_fn : <path> file containing genomic intervals the GATK GenotypeGVCFs command operates over in the genotype_vcf_at_given_alleles.nf process.

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

  // check parameters provided
  validate_parameters()

  // -- MAIN-EXECUTION ------------------------------------------------------
  // prepare panel resource channels 
  PARSE_PANEL_SETTINGS(params.panels_settings)

  reference_ch = PARSE_PANEL_SETTINGS.out.reference_ch // tuple([fasta], panel_name, [fasta_idx])
  annotations_ch = PARSE_PANEL_SETTINGS.out.annotations_ch

  if (params.execution_mode == "in-country") {
    // process in country entry point
    IN_COUNTRY(reference_ch)
    bam_files_ch = IN_COUNTRY.out.bam_files_ch
    sample_tag_reference_files_ch = IN_COUNTRY.out.sample_tag_reference_files_ch
  }

  if (params.execution_mode == "irods") {
    // process IRODS entry point
    IRODS(params.irods_manifest, reference_ch)
    // setup channels for downstream processing
    bam_files_ch = IRODS.out.bam_files_ch // tuple (sample_tag, bam_file, run_id)
    sample_tag_reference_files_ch = IRODS.out.sample_tag_reference_files_ch
  }

  if (params.execution_mode == "aligned_bams"){
    // get bam files channel
    mnf_ch = Channel.fromPath(params.aligned_bams_mnf, checkIfExists: true)
                        | splitCsv(header:true, sep:',')
    bam_files_ch = mnf_ch | map {row -> tuple(row.sample_tag, row.bam_file, row.bam_idx)}
    
    // get sample to pannels relationship channel
    ref_to_sample = mnf_ch | map {row -> tuple(row.panel_name, row.sample_tag)} 

    reference_ch
      | map {it -> tuple(it[1],it[0][0], it[2])} // tuple(panel_name, fasta_file, [fasta_idxs])
      | join(ref_to_sample) // tuple(panel_name, fasta_file, [fasta_idxs], sample_tag)
      | map {it -> tuple(it[3],it[1],it[2],it[0])} // tuple(sample_tag,)
      | set {sample_tag_reference_files_ch}
  }

  COMMON(bam_files_ch, sample_tag_reference_files_ch, annotations_ch)

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
