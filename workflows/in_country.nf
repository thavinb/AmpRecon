#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows
include { bcl_to_cram } from './pipeline-subworkflows/bcl-to-cram.nf'
include { cram_to_bam } from './pipeline-subworkflows/cram-to-bam.nf'

// - process to extract and validate information expected based on input params
include { validate_parameters; load_steps_to_run } from './pipeline-subworkflows/inputHandling.nf'
include { get_taglist_file } from '../modules/manifest2tag.nf'
include { make_samplesheet_manifest } from '../modules/make_samplesheet_manifest.nf'
include { validate_samplesheet_manifest } from '../modules/samplesheet_manifest_validation.nf'
include { PARSE_PANNEL_SETTINGS } from './parse_pannels_settings.nf'


workflow IN_COUNTRY {

	main:



	  // -- Pre Processing ------------------------------------------------------
	  // validate input
 	  validate_parameters()
  
	/*
  	  reference_ch = Channel.from(
                  [file("${params.reference_dir}/grc1/*.fasta"), "PFA_GRC1_v1.0" , file("${params.reference_dir}/grc1/*.fasta.*")],
                  [file("${params.reference_dir}/grc2/*.fasta"), "PFA_GRC2_v1.0", file("${params.reference_dir}/grc2/*.fasta.*")],
                  [file("${params.reference_dir}/spec/*.fasta"), "PFA_Spec", file("${params.reference_dir}/spec/*.fasta.*")]
                  )
  */

      PARSE_PANNEL_SETTINGS()
      reference_ch = PARSE_PANNEL_SETTINGS.out

      // create input channel
      input_csv_ch = Channel.of(tuple(params.run_id,
                                 params.bcl_dir,
                                 params.lane,
                                 params.study_name,
                                 params.read_group,
                                 params.library))
    

       // validate MiSeq run files and directory structure
    
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
       bcl_to_cram(step1_Input_ch)
       manifest_step1_1_Out_ch = bcl_to_cram.out.multiMap { it -> run_id: it[0]
                                                                  mnf: it[1]}
 	
    // get the relevant sample data from the manifest
       ref_tag = Channel.fromPath("${params.results_dir}/*_manifest.csv").splitCsv(header: ["lims_id", "sims_id", "index", "ref", "barcode_sequence", "well", "plate"], skip: 18).map    { row -> tuple(row.lims_id, row.ref, row.index) }

    // group reference files
        reference_ch = Channel.from(
                    [file("$projectDir/references/grc1/Pf_GRC1v1.0.fasta"), "PFA_GRC1_v1.0", file("$projectDir/references/grc1/Pf_GRC1v1.0.fasta.*")],
                    [file("$projectDir/references/grc2/Pf_GRC2v1.0.fasta"), "PFA_GRC2_v1.0", file("$projectDir/references/grc2/Pf_GRC2v1.0.fasta.*")],
                    [file("$projectDir/references/spec/Spec_v1.0.fasta"), "PFA_Spec", file("$projectDir/references/spec/Spec_v1.0.fasta.*")])

        // assign each sample tag the appropriate set of reference files -> tuple('lims_id#index_', 'path/to/reference/genome, 'path/to/reference/index/files')
        ref_tag.combine(reference_ch,  by: 1).map{it -> tuple(it[1]+"#${it[2]}_", it[3], it[4])}.set{sample_tag_reference_files_ch}
 
  

        csv_ch = manifest_step1_1_Out_ch.mnf
    
        // Stage 1 - Step 2: CRAM to BAM
        cram_to_bam(csv_ch, sample_tag_reference_files_ch)
        step1_2_Out_ch = cram_to_bam.out.bam_ch.multiMap { it -> sample_tag: it[0]
                                                                bam_file: it[1]
                                                               run_id:it[2]  }


	step1_3_In_ch = step1_2_Out_ch


	realignment(step1_3_In_ch.sample_tag,
                        step1_3_In_ch.bam_file,
                        step1_3_In_ch.run_id,
                        sample_tag_reference_files_ch,
                        )




}
