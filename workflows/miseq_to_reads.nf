#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

/*
    | MISEQ_TO_READS |-----------------------------------------
    
    This workflow takes a manifest file and path to a Miseq run 
    directory. Basecalls and intensities from this Illumina Miseq
    run are converted to BAM format and demultiplexed. Adapter 
    sequences are found and moved, alignments undone and the
    reads are output in FASTQ format.
    
    ENA submission ready CRAM files are also produced by this 
    workflow. Reads in the FASTQ files are aligned to a reads to a
    reference genome. The resulting file has its header rewritten,
    is sorted by coordinate and is output in CRAM format.

    A channel containing FASTQ file and associated reference genome is
    is output by the workflow. This is in addition to channels that
    link SNP list, annotation file and sample ID with file ID.
    ------------------------------------------------------------------
*/

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
include { upload_pipeline_output_to_s3 } from '../modules/upload_pipeline_output_to_s3.nf'

workflow MISEQ_TO_READS {
  take:
    manifest_file
    reference_ch // tuple (fasta, panel_name, snp_list)

  main:

    // Validate the supplied manifest file
    panel_names_list = reference_ch.map{it -> it[1].toString()}.collect()
    validate_manifest(manifest_file, panel_names_list)

    // convert basecalls
    basecalls_conversion(params.batch_id, params.bcl_dir, params.ena_study_name)

    // Create tag list file
    create_taglist_file(params.ena_study_name, manifest_file)

    // decode multiplexed bam file
    decode_multiplexed_bam(basecalls_conversion.out, create_taglist_file.out)

    // find adapter contamination in bam
    bam_find_adapter(decode_multiplexed_bam.out.decoded_bam_file)

    // split bam by read group into cram
    split_bam_by_readgroup(bam_find_adapter.out)

    // extract index from CRAM file name and add to channel
    split_bam_by_readgroup.out
      | flatten() // /path/to/<runID>_<lane>#<index>.cram
      | map { it -> tuple(it.simpleName.split('#')[1], it) }
      | set {cram_ch} // tuple (index, cram_file)

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

    updated_cram_ch = file_id_ch.map{it -> tuple(it[0], it[3],)} // tuple (file_id, cram_file)
    reference_fasta_ch = file_id_reference_files_ch.map{it -> tuple(it[0], it[2])} // tuple (file_id, reference_fasta_file)

    // --| DEBUG |--------------------------------------------------
    if (!(params.DEBUG_takes_n_bams == null)){
        updated_cram_ch = updated_cram_ch.take(params.DEBUG_takes_n_bams)
    }
    // -------------------------------------------------------------

    updated_cram_ch // tuple(file_id, cram_file)
      // get panel resource files
      | join(reference_fasta_ch) // tuple (file_id, cram_file, ref_fasta)
      | set{ crams_and_reference_files_ch }

    // Create FASTQ and an ENA submission ready CRAM file 
    cram_to_fastq_and_ena_cram(crams_and_reference_files_ch)

    // upload ENA-ready CRAM files to S3 bucket
    if (params.upload_to_s3){
      cram_to_fastq_and_ena_cram.out.cram.set{output_to_s3}
      upload_pipeline_output_to_s3(output_to_s3, "crams")
    }

    // Output channel: fastq_file assigned its appropriate reference fasta
    fastq_files_ch = cram_to_fastq_and_ena_cram.out.fastq // tuple (file_id, fastq_file, reference_fasta_file)

    // Output channel: file_id linked with lims_id
    file_id_ch.map{it -> tuple(it[0], it[2])}.set{file_id_to_sample_id_ch} // tuple (file_id, lims_id)

  emit:
    fastq_files_ch // tuple (file_id, fastq_file, reference_fasta_file)
    file_id_reference_files_ch // tuple (file_id, panel_name, path/to/reference/genome, snp_list)
    file_id_to_sample_id_ch // tuple (file_id, lims_id)
}


workflow {
    // File required for Miseq to Reads input channel
    channel_data = Channel.fromPath(params.channel_data_file, checkIfExists: true)
      .splitCsv(header: true, sep: '\t')

    miseq_to_reads_parameter_check()

    // Miseq to Reads input channels
    manifest  = Channel.fromPath(params.manifest_path, checkIfExists: true)
    reference_ch = channel_data.map { row -> tuple(row.reference_file, row.panel_name, row.snp_list) }

    MISEQ_TO_READS(manifest, reference_ch)
}

def miseq_to_reads_parameter_check(){
    /*
    This functions counts the number of errors on input parameters exclusively used on Incountry workflow
    
    False for any of those conditions counts as an error.
    
    Returns
    ---
    <int> the number of errors found
    */

    def error = 0

    if (params.batch_id == null){
      log.error("A batch_id parameter must be provided for execution mode '${params.execution_mode}'.")
      error += 1
    }

    if (params.bcl_dir == null){
      log.error("A bcl directory must be specified when the 'execution_mode' parameter is set to 'in-country'.")
      error += 1
    }

    if (params.manifest_path){
      manifest = file(params.manifest_path)
      if (!manifest.exists()){
        log.error("${manifest} does not exist.")
        error += 1
      }
    }

    if (params.ena_study_name == null){
      log.error("A study_name parameter must be provided for execution mode '${params.execution_mode}'.")
      error += 1
    }
    if (error > 0) {
      log.error("Parameter errors were found, the pipeline will not run")
      exit 1
    }
}