#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows
include { BCL_TO_CRAM } from './pipeline-subworkflows/bcl-to-cram.nf'
include { DESIGNATE_PANEL_RESOURCES } from './designate_panel_resources.nf'
include { CRAM_TO_BAM } from './pipeline-subworkflows/cram-to-bam.nf'

// - process to extract and validate information expected based on input params
include { validate_parameters } from './input_handling.nf'
include { get_taglist_file } from '../modules/manifest2tag.nf'
include { make_samplesheet_manifest } from '../modules/make_samplesheet_manifest.nf'
include { validate_samplesheet_manifest } from '../modules/samplesheet_manifest_validation.nf'
include { PARSE_PANEL_SETTINGS } from './parse_panels_settings.nf'
include { miseq_run_validation } from '../modules/miseq_run_validation.nf'
include { retrieve_miseq_run_from_s3 } from '../modules/retrieve_miseq_run_from_s3.nf'


workflow IN_COUNTRY {
   take:
      reference_ch // tuple (fasta, panel_name, snp_list)
   
   main:
      
      if (params.s3_bucket_input != null){
         retrieve_miseq_run_from_s3(params.s3_uuid)
         input_csv_ch = retrieve_miseq_run_from_s3.out.tuple_ch
      }
      
      else{
	   //create input channel if not running s3 entrypoint
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
      panel_names_list = reference_ch.map{it -> it[1].toString()}.collect()

      validate_samplesheet_manifest(make_samplesheet_manifest.out.tuple, panel_names_list)

      get_taglist_file_In_ch = input_csv_ch.join(validate_samplesheet_manifest.out)
      get_taglist_file(get_taglist_file_In_ch)

      step1_Input_ch = input_csv_ch.join(get_taglist_file.out)

      // Stage 1 - Step 1: BCL to CRAM
      BCL_TO_CRAM(step1_Input_ch, make_samplesheet_manifest.out.manifest_file)
      cram_ch = BCL_TO_CRAM.out // tuple (sample_tag, cram_fl, run_id)

      // get the relevant sample data from the manifest
      sample_tag_ch = make_samplesheet_manifest.out.manifest_file // WARN: this need to be removed, we should no rely on results dir
                  | splitCsv(header: ["lims_id", "sims_id", "index", "assay",
                                     "barcode_sequence", "well", "plate"],
                                    skip: 18)
                  | map { row -> tuple(row.lims_id, row.assay, row.index) }
                  | map{it -> tuple("${params.run_id}_${params.lane}#${it[2]}_${it[0]}", it[1], it[0])} // tuple (sample_tag, panel_name, lims_id)

      // link file_id to sample_id
      sample_tag_ch.map{it -> tuple("${it[0]}_${it[1]}", it[2])}.set{file_id_to_sample_id_ch}

      // add panel names to sample_tags
      panel_name_cram_ch = cram_ch.join(sample_tag_ch) // tuple (sample_tag, cram_fl, run_id, panel_name)
                           | map{it -> tuple("${it[0]}_${it[3]}", it[1], it[2])} // tuple(sample_tag_panel_name, cram_fl, run_id)
      new_sample_tag_ch = sample_tag_ch.map{it -> tuple("${it[0]}_${it[1]}", it[1])} // tuple (new_sample_tag, panel_name)  

      // assign each sample tag the appropriate set of reference files
      DESIGNATE_PANEL_RESOURCES(new_sample_tag_ch, reference_ch)
      sample_tag_reference_files_ch = DESIGNATE_PANEL_RESOURCES.out.sample_tag_reference_files_ch 
      // tuple(new_sample_tag, panel_name path/to/reference/genome, snp_list)

      // Stage 1 - Step 2: CRAM to BAM
      CRAM_TO_BAM(panel_name_cram_ch, sample_tag_reference_files_ch.map{it -> tuple(it[0], it[2], it[1])})  // tuple (new_sample_tag, ref_fasta, panel_name)
      bam_files_ch = CRAM_TO_BAM.out.bam_ch

   emit:
      bam_files_ch // tuple (new_sample_tag, bam_file)
      sample_tag_reference_files_ch // tuple (new_sample_tag, panel_name, path/to/reference/genome, snp_list)
      file_id_to_sample_id_ch // tuple (file_id, sample_id)
}
