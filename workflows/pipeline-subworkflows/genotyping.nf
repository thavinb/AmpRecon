#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { gatk_haplotype_caller_gatk4 } from '../../modules/gatk_haplotype_caller_gatk4.nf'
include { bqsr } from '../../modules/bqsr.nf' addParams(gatk:params.gatk3)
include { genotype_vcf_at_given_alleles } from '../../modules/genotype_vcf_at_given_alleles.nf' addParams(gatk:params.gatk3, bgzip:'bgzip')

workflow GENOTYPING {
  //

  take:
        input_sample_tags_bams_indexes
        sample_tag_reference_files_ch
  main:

    // base quality score recalibration
    input_sample_tags_bams_indexes.map{it -> it}.join(sample_tag_reference_files_ch).map{ it ->
        tuple(it[0], it[1], it[2], it[3], it[3]+".fai", it[3]+".dict")}.set{bqsr_input}
    bqsr(bqsr_input)

    // haplotype caller
    bqsr.out.join(sample_tag_reference_files_ch).map{ it ->
    tuple(it[0], it[1], it[2], it[2]+".fai", it[2]+".dict")}.set{haplotype_caller_input}
    gatk_haplotype_caller_gatk4(haplotype_caller_input)

    // genotype alleles in VCFs
    gatk_haplotype_caller_gatk4.out.vcf_file_and_index.join(sample_tag_reference_files_ch).map{ it ->
        tuple(it[0],it[1], it[2], it[3], it[3]+".fai", it[3]+".dict")}.set{genotyping_input_ch}   
    genotype_vcf_at_given_alleles(genotyping_input_ch)

//  emit:
}
