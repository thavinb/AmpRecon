#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bcftools_mpileup } from '../../modules/bcftools_genotyping.nf'
include { bcftools_call } from '../../modules/bcftools_genotyping.nf'
include { bcftools_filter } from '../../modules/bcftools_genotyping.nf'
include { upload_pipeline_output_to_s3 } from '../../modules/upload_pipeline_output_to_s3.nf'

workflow GENOTYPING_BCFTOOLS {

  take:
        input_sample_tags_bams_indexes // tuple(sample_tag, bam_file, bam_index_file)
        sample_tag_reference_files_ch // tuple(sample_tag, reference_fasta, reference_fasta_index, reference_ploidy_file, reference_annotation_file)
  main:

    // compute genotype likelihoods
    input_sample_tags_bams_indexes // tuple (sample_tag, bam_file, bam_index)
        | join(sample_tag_reference_files_ch) // tuple (sample_tag, bam_file, bam_index, reference_fasta, reference_fasta_index, reference_ploidy_file, reference_annotation_file)
        | map{ it -> tuple(it[0], it[1], it[2], it[3], it[4], it[6])}
        | set{mpileup_input} // tuple(sample_tag, bam_file, bam_index, reference_fasta, reference_fasta_index, reference_annotation_vcf)
    bcftools_mpileup(mpileup_input)

    // call SNP sites
    bcftools_mpileup.out // tuple (sample_tag, bcf_file)
	| join(sample_tag_reference_files_ch) // tuple (sample_tag, bcf_file, reference_fasta, reference_fasta_index, reference_ploidy_file, reference_annotation_file)
	| map{ it -> tuple(it[0], it[1], it[4])}
	| set{call_input} // tuple (sample_tag, bcf_file, reference_ploidy_file)
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
