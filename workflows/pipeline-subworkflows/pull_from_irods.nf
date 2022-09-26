#!/usr/bin/env nextflow
// enable dsl2
nextflow.enable.dsl = 2

// import irods processes
include { irods_manifest_parser } from '../../modules/irods_manifest_parser.nf'
include { irods_retrieve } from '../../modules/irods_retrieve.nf'
include { scramble_cram_to_bam } from '../../modules/scramble.nf'


workflow PULL_FROM_IRODS {
  
  take:
    irods_ch // tuple(id_run, WG_lane)
    sample_id_ref_ch // tuple(WG_lane, run_id, fasta_file, fasta_idx_files, pannel_name) * needed to associate new_ids to pannels
  
  main:
    // get names and paths 
    irods_manifest_parser(irods_ch) // tuple(irods_sample_id, irods_flpath, WG_lane)
    

    // get new_sample_ids relationship with pannels resources
    irods_manifest_parser.out // tuple(new_sample_id, iRODS_file_path, WG_lane)
              | map {it -> tuple( it[2], "${it[2]}_${it[0]}")}
=======

    // get new_sample_ids relationship with pannels resources
    irods_manifest_parser.out // tuple(new_sample_id, iRODS_file_path, id_run, WG_lane)
              | map {it -> tuple( it[3], "${it[3]}_${it[0]}-")}
              | set {newSample_WgLn_ch} // tuple(WG_lane, new_sample_id)
    
    sample_id_ref_ch
              // no need for run_id on this channel
              | map { it -> tuple( it[0], it[2], it[3], it[4])} // tuple(WG_lane, fasta_file, fasta_idx_files)
              // associate with new_sample_ids
              | join (newSample_WgLn_ch)  // tuple(WG_lane, fasta_file, fasta_idx_files, pannel_name, new_sample_id)
              | set {old_new_and_pannel_ch}
    
    old_new_and_pannel_ch
              | map { it -> tuple("${it[4]}${it[3]}", it[1], it[2], it[3],)} // tuple( new_sample_id_pannel, fasta_file, fasta_idx_files, pannel_name)
              | set { sample_tag_reference_files_ch }
    
    // add new_sample_id to irods_retrieve
    //add new id key [WGlane]_[sample_id]-[pannel]
    irods_manifest_parser.out
        | map {it -> tuple(it[3], "${it[3]}_${it[0]}-", it[1], it[2])} // tuple (WG_lane, [WG_lane]_[irods_sample_id]-, irods_flpath, id_run)
        | join (old_new_and_pannel_ch) // tuple (WG_lane, new_sample_id, irods_flpath, id_run, fasta_file, fasta_idx_files, pannel_name, new_sample_id)
        | map {it -> tuple("${it[1]}${it[6]}", it[2], it[3])}
        | set {irods_retrieve_In_ch} // tuple (new_sample_id_pannel, irods_flpath, id_run)

    sample_tag_reference_files_ch.first().view()
    // Retrieve CRAM files from iRODS
    irods_retrieve(irods_retrieve_In_ch)

    // Convert iRODS CRAM files to BAM format
    scramble_cram_to_bam(irods_retrieve.out)
    bam_files_ch = scramble_cram_to_bam.out

    // --------------------------------------------------------------------

  emit:
    bam_files_ch  // tuple(new_sample_id, bam_file)
    sample_tag_reference_files_ch // tuple(new_sample_id, ref_files)
}

/*
// -------------------------- DOCUMENTATION -----------------------------------
*/
