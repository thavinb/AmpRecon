#!/usr/bin/env nextflow
// enable dsl2
nextflow.enable.dsl = 2

// import irods processes
include { irods_retrieve } from '../../modules/irods_retrieve.nf'
include { scramble_cram_to_bam } from '../../modules/scramble.nf'


workflow PULL_FROM_IRODS {
  
  take:
    irods_ch // tuple(new_sample_id, irods_path)
  
  main:
    // Retrieve CRAM files from iRODS
    irods_retrieve(irods_ch)

    // Convert iRODS CRAM files to BAM format
    scramble_cram_to_bam(irods_retrieve.out)
    bam_files_ch = scramble_cram_to_bam.out

    // --------------------------------------------------------------------
  emit:
    bam_files_ch  // tuple(new_sample_id, bam_file)
}

/*
// -------------------------- DOCUMENTATION -----------------------------------
*/
