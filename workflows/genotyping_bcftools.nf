#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bcftools_mpileup } from '../modules/bcftools_genotyping.nf'
include { bcftools_call } from '../modules/bcftools_genotyping.nf'
include { bcftools_filter } from '../modules/bcftools_genotyping.nf'
include { upload_pipeline_output_to_s3 } from '../modules/upload_pipeline_output_to_s3.nf'

workflow GENOTYPING_BCFTOOLS {
  take:
    input_file_ids_bams_indexes // tuple(file_id, bam_file, bam_index_file)
    file_id_reference_files_ch // tuple(file_id, reference_fasta, snp_list)

  main:

    // compute genotype likelihoods
    input_file_ids_bams_indexes // tuple (file_id, bam_file, bam_index)
      | join(file_id_reference_files_ch) // tuple (file_id, bam_file, bam_index, reference_fasta, snp_list)
      | set{mpileup_input} // tuple(file_id, bam_file, bam_index, reference_fasta, snp_list)
    bcftools_mpileup(mpileup_input)

    // call SNP sites
    bcftools_mpileup.out // tuple (file_id, bcf_file)
	| join(file_id_reference_files_ch) // tuple (file_id, bcf_file, reference_file, snp_list)
	| map{ it -> tuple(it[0], it[1])}
	| set{call_input} // tuple (file_id, bcf_file)
    bcftools_call(call_input)

    // filter genetic variants
    bcftools_filter(bcftools_call.out)
    bcftools_filter.out.set{genotyped_vcf_ch}

    // upload VCF files to S3 bucket
    if (params.upload_to_s3){
      bcftools_filter.out.map{ it -> it[1] }.set{output_to_s3}
      upload_pipeline_output_to_s3(output_to_s3)
    }

  emit:
    genotyped_vcf_ch
}

