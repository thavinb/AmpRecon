#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { gatk_haplotype_caller_gatk4 } from '../../modules/gatk_haplotype_caller_gatk4.nf'
include { genotype_vcf_at_given_alleles } from '../../modules/genotype_vcf_at_given_alleles.nf' addParams(gatk:params.gatk3, bgzip:'bgzip')
include { index_gzipped_vcf } from '../../modules/index_gzipped_vcf.nf'
include { upload_pipeline_output_to_s3 } from '../../modules/upload_pipeline_output_to_s3.nf'

workflow GENOTYPING_GATK {

  take:
        input_sample_tags_bams_indexes
        sample_tag_reference_files_ch
  main:

    input_sample_tags_bams_indexes // tuple (sample_tag, bam_file, bam_index)
            | join(sample_tag_reference_files_ch) // tuple (sample_tag, bam_file, bam_index [ref_files])
            | map{ it -> tuple(it[0], it[1], it[2], it[3], it[3]+".fai", it[3]+".dict")}
            | set{haplotype_caller_input} // tuple(sample_tag, bam_file, bam_index, [ref_files], )

    gatk_haplotype_caller_gatk4(haplotype_caller_input)

    // genotype alleles in VCFs
    gatk_haplotype_caller_gatk4.out.vcf_file_and_index
        | join(sample_tag_reference_files_ch)
        | map{ it -> tuple(it[0],it[1], it[2], it[3], it[3]+".fai", it[3]+".dict")}
        | set{genotyping_input_ch}   
    
    genotype_vcf_at_given_alleles(genotyping_input_ch)
    index_gzipped_vcf(genotype_vcf_at_given_alleles.out).set{genotyped_vcf_ch}

    // upload VCF files / indices to S3 bucket
    if (params.upload_to_s3){
      index_gzipped_vcf.out.map{it -> tuple(it[1], it[2])}.flatten().set{output_to_s3}
      upload_pipeline_output_to_s3(output_to_s3)
    }    

  emit:
    genotyped_vcf_ch
}
