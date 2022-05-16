#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include { bcl_to_cram } from '../bespoke_modules/bcl_to_cram.nf'
include { fetch_from_irods } from '../bespoke_modules/fetch_from_irods.nf'
include { cram_to_bam } from '../bespoke_modules/cram_to_bam.nf'

workflow reads_retrieval_workflow {
    take:
        bcl_dir
        irods_paths
    main:
        bcl_to_cram (bcl_dir)
        bcl_to_cram_ch = bcl_to_cram.out

        fetch_from_irods (irods_paths)
        irods_ch = fetch_from_irods.out

        bcl_to_cram_ch.concat(irods_ch).set{cram_ch}

        cram_to_bam (cram_ch)
        bam_ch = cram_to_bam.out
    emit:
        bam_ch
}


