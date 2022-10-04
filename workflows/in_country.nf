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
      reference_ch // tuple ([fasta_file], pannel_name, [fasta_idxs])
   
   main:
      
      if (params.uuid != null){
         retrieve_miseq_run_from_s3(params.uuid)
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
      validate_samplesheet_manifest(make_samplesheet_manifest.out.tuple)

      get_taglist_file_In_ch = input_csv_ch.join(validate_samplesheet_manifest.out)
      get_taglist_file(get_taglist_file_In_ch)

      step1_Input_ch = input_csv_ch.join(get_taglist_file.out)

      // Stage 1 - Step 1: BCL to CRAM
      BCL_TO_CRAM(step1_Input_ch, make_samplesheet_manifest.out.manifest_file)
      cram_ch = BCL_TO_CRAM.out // tuple (sample_tag, cram_fl, run_id)
      
      //manifest_step1_1_Out_ch = BCL_TO_CRAM.out.multiMap { it -> run_id: it[0]
      //                                                           mnf: it[1]}

      // get the relevant sample data from the manifest
      ref_tag = make_samplesheet_manifest.out.manifest_file // WARN: this need to be removed, we should no rely on results dir
                  | splitCsv(header: ["lims_id", "sims_id", "index", "assay",
                                     "barcode_sequence", "well", "plate"],
                                    skip: 18)
                  | map { row -> tuple(row.lims_id, row.assay, row.index) }

      // assign each sample tag the appropriate set of reference files -> tuple('lims_id#index_', 'path/to/reference/genome, 'path/to/reference/index/files')
      
      ref_tag // tuple (lims_id, pannel_name, index)
         | combine(reference_ch,  by: 1) // tuple (pannel_name, lims_id, index, [fasta_file], [fasta_idxs])
         | map{it -> tuple("${params.run_id}_${params.lane}#${it[2]}_${it[1]}", it[3][0], it[4], it[0])}
         | set{sample_tag_reference_files_ch}

      //csv_ch = manifest_step1_1_Out_ch.mnf

      // Stage 1 - Step 2: CRAM to BAM
      CRAM_TO_BAM(cram_ch, sample_tag_reference_files_ch)
      bam_files_ch = CRAM_TO_BAM.out.bam_ch

   emit:
      bam_files_ch // tuple (sample_tag, bam_file)
      sample_tag_reference_files_ch // tuple('lims_id#index_', 'path/to/reference/genome, 'path/to/reference/index/files', pannel_name)
}
