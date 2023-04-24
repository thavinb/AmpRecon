#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include { basecalls_conversion } from '../modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from '../modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'
include { rename_cram_fls } from '../modules/rename_cram_fls.nf'
include { cram_to_fastq_and_ena_cram } from '../modules/cram_to_fastq_and_ena_cram.nf'

// - process to extract and validate information expected based on input params
include { get_taglist_file } from '../modules/manifest2tag.nf'
include { make_samplesheet_manifest } from '../modules/make_samplesheet_manifest.nf'
include { validate_samplesheet_manifest } from '../modules/samplesheet_manifest_validation.nf'
include { retrieve_miseq_run_from_s3 } from '../modules/retrieve_miseq_run_from_s3.nf'

workflow BCL_TO_COLLATED_CRAM {
    take:
        run_id
        bcl_dir
        lane
        study_name
        read_group
        library
        tag_file
        manifest

    main:
        // convert basecalls
        basecalls_conversion(run_id, bcl_dir, lane, study_name, read_group)

        // decode multiplexed bam file
        decode_multiplexed_bam(basecalls_conversion.out, tag_file)

        // find adapter contamination in bam
        bam_find_adapter(decode_multiplexed_bam.out.decoded_bam_file)
  
        // split bam by read group into cram
        bam_to_cram(bam_find_adapter.out)

        // rename samples to samplesheet provided names
        // on this step the pattern for file names are set as
        // [run_id]_[lane]#[index]_[sample_name]
        rename_cram_fls(run_id,
                        decode_multiplexed_bam.out.bam_metrics_file,
                        bam_to_cram.out,
                        lane,
                        manifest
                        )
        cram_ch = rename_cram_fls.out

        // add new sample tag to cram_ch
        cram_ch
          | flatten()
          | map { it -> tuple(it.simpleName, it) }
          | set {final_cram_ch} // tuple (file_id, cram_file,)

    emit:
        final_cram_ch // tuple (file_id, cram_file)

}

workflow COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ {
/*
    | COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ |-----------------------------------------
    
    This workflow does the following to each CRAM file it is supplied:
    [1] Undoes alignments, moves the adapter sequences and outputs the reads in fastq format.
    [2] Aligns reads to a reference, rewrites the header, sorts by coordinate and outputs the alignments in CRAM format.
    ------------------------------------------------------------------
*/
    take:
        cram_ch // tuple (file_id, cram_file)
        file_id_reference_files_ch // tuple (file_id, ref_fasta)

    main:
        // --| DEBUG |--------------------------------------------------
        if (!(params.DEBUG_takes_n_bams == null)){
            cram_ch = cram_ch.take(params.DEBUG_takes_n_bams)
        }
        // -------------------------------------------------------------

        cram_ch.first().view()
        file_id_reference_files_ch.first().view()
        cram_ch // tuple(file_id, cram_file)
              // get panel resource files
              | join(file_id_reference_files_ch) // tuple (file_id, cram_file, ref_fasta)
              | set{ crams_and_reference_files_ch }

        // Create FASTQ and an ENA submission ready CRAM file 
        cram_to_fastq_and_ena_cram(crams_and_reference_files_ch)
        fastq_ch = cram_to_fastq_and_ena_cram.out.fastq
  
    emit:
        fastq_ch  // tuple (file_id, fastq_file, reference_fasta_file, panel_name)
}

workflow MISEQ_TO_READS {
  take:
    reference_ch // tuple (fasta, panel_name, snp_list)

  main:
      
    if (params.s3_bucket_input != null) {
      retrieve_miseq_run_from_s3(params.s3_uuid)
      bcl_dir = retrieve_miseq_run_from_s3.out
    }
  
    else {
	    // Use supplied bcl dir if not running s3 entrypoint
      bcl_dir = params.bcl_dir
    }
  
    // process samplesheets manifest (necessary to get barcodes) and validate it
    make_samplesheet_manifest(params.run_id, params.bcl_dir)
    manifest_file = make_samplesheet_manifest.out
    panel_names_list = reference_ch.map{it -> it[1].toString()}.collect()

    manifest_file.view()
    panel_names_list.view()
    // Validate the manifest created from the samplesheet file
    validate_samplesheet_manifest(manifest_file, panel_names_list)

    //
    get_taglist_file(params.run_id,
                    params.study_name,
                    params.library,
                    manifest_file)

    // Stage 1 - Step 1: BCL to CRAM
    BCL_TO_COLLATED_CRAM(params.run_id,
                        bcl_dir,
                        params.lane, 
                        params.study_name, 
                        params.read_group, 
                        params.library, 
                        get_taglist_file.out, 
                        manifest_file)
    cram_ch = BCL_TO_COLLATED_CRAM.out // tuple (file_id, cram_file, run_id)

    // get the relevant sample data from the manifest
    file_id_ch = manifest_file
      | splitCsv(header: ["lims_id", "sims_id", "index", "assay",
          "barcode_sequence", "well", "plate"],
          skip: 18)
      | map { row -> tuple(row.lims_id, row.assay, row.index) }
      | map{it -> tuple("${params.run_id}_${params.lane}#${it[2]}_${it[0]}", it[1], it[0])} // tuple (file_id, panel_name, lims_id)

    // link file_id to sample_id
    file_id_ch.map{it -> tuple("${it[0]}_${it[1]}", it[2])}.set{file_id_to_sample_id_ch}

    // add panel names to file_ids
    panel_name_cram_ch = cram_ch.join(file_id_ch) // tuple (file_id, cram_file, panel_name)
     | map{it -> tuple("${it[0]}_${it[2]}", it[1])} // tuple (file_id_panel_name, cram_file)
    new_file_id_ch = file_id_ch.map{it -> tuple("${it[0]}_${it[1]}", it[1])} // tuple (file_id, panel_name)  

    // assign each sample tag the appropriate set of reference files
    new_file_id_ch
      | combine(reference_ch,  by: 1) // tuple (panel_name, file_id, fasta,snp_list)
      | map{it -> tuple(it[1], it[0], it[2], it[3])}
      | set{file_id_reference_files_ch}
    // tuple(file_id, panel_name path/to/reference/genome, snp_list)

    // Stage 1 - Step 2: Collated CRAM to FASTQ and split CRAM
    COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ(panel_name_cram_ch, file_id_reference_files_ch.map{it -> tuple(it[0], it[2])})  // tuple (new_file_id, ref_fasta)
    fastq_files_ch = COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ.out.fastq_ch

  emit:
    fastq_files_ch // tuple (file_id, fastq_file, reference_fasta_file)
    file_id_reference_files_ch // tuple (file_id, panel_name, path/to/reference/genome, snp_list)
    file_id_to_sample_id_ch // tuple (file_id, sample_id)
}


def miseq_to_reads_parameter_check(){
    /*
    This functions counts the number of errors on input parameters exclusively used on Incountry workflow
    
    checks:
     - if irods manifest was provided
     - if irods manifest provided exists
    
    False for any of those conditions counts as an error.
    
    Returns
    ---
    <int> the number of errors found
    */

    def err = 0

    if (params.run_id == null){
      log.error("A run_id parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }

    if (params.bcl_dir == null && params.s3_bucket_input == null){
      log.error("Either a bcl directory or a s3 bucket input must be specified for in-country execution_mode.")
      err += 1
    } 
    

    if (params.lane == null){
      log.error("A lane parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    
    if (params.study_name == null){
      log.error("A study_name parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    if (params.read_group == null){
      log.error("A read_group parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    
    if (params.library == null){
      log.error("A library parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    return err
}

