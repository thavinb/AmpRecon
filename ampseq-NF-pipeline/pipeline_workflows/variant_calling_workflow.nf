#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { variant_calling  } from '../bespoke_modules/variant_calling.nf'

workflow variant_calling_workflow {
    take:
        filtered_bam_ch
    main:
        variant_calling (filtered_bam_ch)
        variant_calling.out.set { vcf_ch}
    emit:
        vcf_ch
}
