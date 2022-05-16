#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include { bcl_to_bam } from '../bespoke_modules/bcl_to_bam.nf'
include { fetch_from_irods } from '../bespoke_modules/fetch_from_irods.nf'

workflow reads_retrieval_workflow {
    take:
        bcl_dir
        irods_paths
    main:
        bcl_to_bam (bcl_dir)
        bcl_to_bam_ch = bcl_to_bam.out
        fetch_from_irods (irods_paths)
        irods_ch = fetch_from_irods.out
        bcl_to_bam_ch.join(irods_ch).set{bam_ch}
    emit:
        bam_ch
}


