#!/usr/bin/env nextflow
// enable dsl2
nextflow.enable.dsl = 2

// import irods processes
//include { irods_manifest_parser } from '../../modules/irods_manifest_parser.nf'
include { irods_retrieve } from '../../modules/irods_retrieve.nf'
include { scramble_cram_to_bam } from '../../modules/scramble.nf'


workflow PULL_FROM_IRODS {
  
  take:
    irods_ch // tuple(sample_id, WG_lane, irods_path)
    sample_id_ref_ch // tuple(WG_lane, run_id, fasta_file, fasta_idx_files, panel_name) * needed to associate new_ids to panels
  
  main:
    // get new_sample_ids relationship with panels resources
    irods_ch // tuple(sample_id, WG_lane, irods_path)
        | map {it -> tuple( it[1], "${it[1]}_${it[0]}_")}
        | set {newSample_WgLn_ch} // tuple(WG_lane, new_sample_id)

    sample_id_ref_ch
              // no need for run_id on this channel
              | map { it -> tuple( it[0], it[2], it[3], it[4])} // tuple(WG_lane, fasta_file, fasta_idx_files)
              // associate with new_sample_ids
              | join (newSample_WgLn_ch)  // tuple(WG_lane, fasta_file, fasta_idx_files, panel_name, new_sample_id)
              | set {old_new_and_panel_ch}

    old_new_and_panel_ch
              | map { it -> tuple("${it[4]}${it[3]}", it[1], it[2], it[3],)} // tuple( new_sample_id_panel, fasta_file, fasta_idx_files, panel_name)
              | set { sample_tag_reference_files_ch }
    
    // add new_sample_id to irods_retrieve
    // add new id key [WGlane]_[sample_id]-[panel]
    irods_ch
        | map {it -> tuple(it[1], "${it[1]}_${it[0]}_", it[2])}//, it[2])} // tuple (WG_lane, [WG_lane]_[irods_sample_id]-, irods_flpath)
        | join (old_new_and_panel_ch) // tuple (WG_lane, new_sample_id, irods_flpath, fasta_file, fasta_idx_files, panel_name, new_sample_id)
        | map {it -> tuple("${it[1]}${it[5]}", it[2])}//, it[3])}
        | set {irods_retrieve_In_ch} // tuple (new_sample_id_panel, irods_flpath)
    
    // Retrieve CRAM files from iRODS
    irods_retrieve(irods_retrieve_In_ch)
    // Convert iRODS CRAM files to BAM format
    scramble_cram_to_bam(irods_retrieve.out)
    bam_files_ch = scramble_cram_to_bam.out

    // --------------------------------------------------------------------
  emit:
    bam_files_ch  // tuple(new_sample_id, bam_file)
    sample_tag_reference_files_ch // tuple(new_sample_id, fasta_file, ref_files, fasta_index, panel_names )
}

/*
// -------------------------- DOCUMENTATION -----------------------------------
*/
