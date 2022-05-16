#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include { variant_calling  } from '../bespoke_modules/variant_calling.nf'

workflow variant_calling_workflow {
    take:
        bam_ch
    main:
        variant_calling (bam_ch)
    emit:
        variant_calling.out
}
