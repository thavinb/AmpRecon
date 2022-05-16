#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { genotype_amino_alleles  } from '../bespoke_modules/genotype_amino_alleles.nf'
include { assemble_amino_haplotypes  } from '../bespoke_modules/assemble_amino_haplotypes.nf'

workflow genotyping_workflow {
    take:
        vcf_ch
    main:
        genotype_amino_alleles (vcf_ch)
        assemble_amino_haplotypes (genotype_amino_alleles.out)
    emit:
        genotype_amino_alleles.out
        assemble_amino_haplotypes.out
}
