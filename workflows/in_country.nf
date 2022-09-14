#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows
include { BCL_TO_CRAM } from './pipeline-subworkflows/bcl-to-cram.nf'
include { CRAM_TO_BAM } from './pipeline-subworkflows/cram-to-bam.nf'

// - process to extract and validate information expected based on input params
include { validate_parameters } from './input_handling.nf'
include { get_taglist_file } from '../modules/manifest2tag.nf'
include { make_samplesheet_manifest } from '../modules/make_samplesheet_manifest.nf'
include { validate_samplesheet_manifest } from '../modules/samplesheet_manifest_validation.nf'
include { PARSE_PANNEL_SETTINGS } from './parse_pannels_settings.nf'
include { miseq_run_validation } from '../modules/miseq_run_validation.nf'
include { retrieve_miseq_run_from_s3 } from '../modules/retrieve_miseq_run_from_s3.nf'


workflow IN_COUNTRY {
   take:
      reference_ch
   
   main:
      if ( params.download_from_s3 == true & !file("${params.bcl_dir}").exists() ) {
         retrieve_miseq_run_from_s3(params.bcl_id)
         input_csv_ch = retrieve_miseq_run_from_s3.out.tuple_ch
      } else {
         // create input channel
         input_csv_ch = Channel.of(tuple(params.run_id,
                                 params.bcl_dir,
                                 params.lane,
                                 params.study_name,
                                 params.read_group,
                                 params.library))
      }

      // process samplesheets manifest (necessary to get barcodes) and validate it
      input_csv_ch
         | map {it -> tuple (it[0], it[1])} // tuple(run_id, bcl_dir)
         | set { make_samplesheet_In_ch}

      make_samplesheet_manifest(make_samplesheet_In_ch)
      validate_samplesheet_manifest(make_samplesheet_manifest.out.tuple)

      get_taglist_file_In_ch = input_csv_ch.join(validate_samplesheet_manifest.out)
      get_taglist_file(get_taglist_file_In_ch)

      step1_Input_ch = input_csv_ch.join(get_taglist_file.out)

      // Stage 1 - Step 1: BCL to CRAM
      BCL_TO_CRAM(step1_Input_ch)
      manifest_step1_1_Out_ch = BCL_TO_CRAM.out.multiMap { it -> run_id: it[0]
                                                                 mnf: it[1]}
 	
      // get the relevant sample data from the manifest
      ref_tag = Channel.fromPath("${params.results_dir}/*_manifest.csv")
                  | splitCsv(header: ["lims_id", "sims_id", "index", "ref",
                                     "barcode_sequence", "well", "plate"],
                                    skip: 18)
                  | map { row -> tuple(row.lims_id, row.ref, row.index) }

      // assign each sample tag the appropriate set of reference files -> tuple('lims_id#index_', 'path/to/reference/genome, 'path/to/reference/index/files')
      ref_tag | combine(reference_ch,  by: 1)
              | map{it -> tuple(it[1]+"#${it[2]}_", it[3], it[4])}
              | set{sample_tag_reference_files_ch}
 
      csv_ch = manifest_step1_1_Out_ch.mnf
    
      // Stage 1 - Step 2: CRAM to BAM
      CRAM_TO_BAM(csv_ch, sample_tag_reference_files_ch)
      bam_files_ch = CRAM_TO_BAM.out.bam_ch 

   emit:
      bam_files_ch // tuple (sample_tag, bam_file)
      sample_tag_reference_files_ch // tuple('lims_id#index_', 'path/to/reference/genome, 'path/to/reference/index/files')
}
