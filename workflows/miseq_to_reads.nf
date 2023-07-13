#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include { basecalls_conversion } from '../modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from '../modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { split_bam_by_readgroup } from '../modules/split_bam_by_readgroup.nf'
include { cram_to_fastq_and_ena_cram } from '../modules/cram_to_fastq_and_ena_cram.nf'

// - process to extract and validate information expected based on input params
include { create_taglist_file } from '../modules/create_taglist_file.nf'
include { validate_manifest } from '../modules/validate_manifest.nf'
include { retrieve_miseq_run_from_s3 } from '../modules/retrieve_miseq_run_from_s3.nf'

workflow BCL_TO_COLLATED_CRAM {
    take:
        run_id
        bcl_dir
        study_name
        manifest_file

    main:
        // convert basecalls
        basecalls_conversion(run_id, bcl_dir, study_name)

        // Create tag list file
        create_taglist_file(study_name, manifest_file)

        // decode multiplexed bam file
        decode_multiplexed_bam(basecalls_conversion.out, create_taglist_file.out)

        // find adapter contamination in bam
        bam_find_adapter(decode_multiplexed_bam.out.decoded_bam_file)
  
        // split bam by read group into cram
        split_bam_by_readgroup(bam_find_adapter.out)
        cram_ch = split_bam_by_readgroup.out

        // extract index from CRAM file name and add to channel
        cram_ch
          | flatten() // /path/to/runID_lane#index.cram
          | map { it -> tuple(it.simpleName.split('#')[1], it) }
          | set {final_cram_ch} // tuple (index, cram_file,)

    emit:
        final_cram_ch // tuple (index, cram_file)

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
    manifest_file
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

    // Validate the supplied manifest file
    panel_names_list = reference_ch.map{it -> it[1].toString()}.collect()
    validate_manifest(manifest_file, panel_names_list)

    // Stage 1 - Step 1: BCL to CRAM
    BCL_TO_COLLATED_CRAM(params.run_id,
                        bcl_dir,
                        params.ena_study_name,
                        manifest_file)
    cram_ch = BCL_TO_COLLATED_CRAM.out // tuple (file_id, cram_file)

    // get the relevant sample data from the manifest and link with cram files
    file_id_ch = manifest_file
      | splitCsv(header: true, sep: '\t')
      | map { row -> tuple(row.barcode_number, row.primer_panel, row.sample_id) }
      // Match (by index) assay and lims_id with associated CRAM file
      | join(cram_ch)  // tuple (index, panel_name, sample_id, cram_file) 
      // append lims_id and assay panel_name to file_id (CRAM file simpleName)
      | map{it -> tuple("${it[3].simpleName}_${it[2]}_${it[1]}", it[1], it[2], it[3])} // tuple (fileId_limsId_panelName, panel_name, lims_id, cram_file)

    // Output channel: file_id assigned appropriate set of assay panel resource files
    file_id_ch
      | map{it -> tuple(it[0], it[1])} // tuple (file_id, panel_name)  
      | combine(reference_ch,  by: 1) // tuple (panel_name, file_id, fasta, snp_list)
      | map{it -> tuple(it[1], it[0], it[2], it[3])}
      | set{file_id_reference_files_ch} // tuple (file_id, panel_name, path/to/reference/genome, snp_list)

    // Stage 1 - Step 2: Collated CRAM to FASTQ and split CRAM
    updated_cram_ch = file_id_ch.map{it -> tuple(it[0], it[3],)} // tuple (file_id, cram_file)
    reference_fasta_ch = file_id_reference_files_ch.map{it -> tuple(it[0], it[2])} // tuple (file_id, reference_fasta_file)
    COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ(updated_cram_ch, reference_fasta_ch)

    // Output channel: fastq_file assigned its appropriate reference fasta
    fastq_files_ch = COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ.out.fastq_ch // tuple (file_id, fastq_file, reference_fasta_file)

    // Output channel: file_id linked with lims_id
    file_id_ch.map{it -> tuple(it[0], it[2])}.set{file_id_to_sample_id_ch} // tuple (file_id, lims_id)

  emit:
    fastq_files_ch // tuple (file_id, fastq_file, reference_fasta_file)
    file_id_reference_files_ch // tuple (file_id, panel_name, path/to/reference/genome, snp_list)
    file_id_to_sample_id_ch // tuple (file_id, lims_id)
}

	
workflow {
    // File required for alignment input channel
    manifest  = Channel.fromPath(params.manifest_path, checkIfExists: true)

    // Input reference channel
    reference_ch = channel_data.map { row -> tuple(row.reference_file, row.panel_name, row.snp_list) }
  
    MISEQ_TO_READS(manifest, reference_ch)
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

    if (params.manifest_path){
      samplesheet = file(params.manifest_path)
      if (!manifest.exists()){
        log.error("${manifest} does not exist.")
      }
    }

    if (params.ena_study_name == null){
      log.error("A study_name parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    return err
}

