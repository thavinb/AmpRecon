#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

// enable dsl2
nextflow.enable.dsl = 2

include { bcftools_mpileup } from '../modules/bcftools_genotyping.nf'
include { bcftools_call } from '../modules/bcftools_genotyping.nf'
include { bcftools_filter } from '../modules/bcftools_genotyping.nf'
include { upload_pipeline_output_to_s3 } from '../modules/upload_pipeline_output_to_s3.nf'

workflow GENOTYPING_BCFTOOLS {
  take:
    input_file_ids_bams_indexes // tuple(file_id, bam_file, bam_index_file)
    file_id_reference_files_ch // tuple(file_id, reference_fasta, snp_list, sample_key)

  main:
     
    // compute genotype likelihoods
    input_file_ids_bams_indexes // tuple (sample_key, bam_file, bam_index)
      | join(file_id_reference_files_ch) // tuple (sample_key, file_id, bam_file, bam_index, reference_fasta, snp_list)
      | map { it -> tuple(it[3], it[1], it[2], it[4], it[5]) }
      | set{mpileup_input} // tuple(file_id, bam_file, bam_index, reference_fasta, snp_list)

    bcftools_mpileup(mpileup_input)

    // call SNP sites 
    bcftools_call(bcftools_mpileup.out)

    // filter genetic variants
    bcftools_filter(bcftools_call.out)
    bcftools_filter.out.set{genotyped_vcf_ch}

    // upload VCF files to S3 bucket
    if (params.upload_to_s3){
      bcftools_filter.out.map{ it -> it[1] }.set{output_to_s3}
      upload_pipeline_output_to_s3(output_to_s3, "vcfs")
    }

  emit:
    genotyped_vcf_ch
}

workflow {
    // File required for BCFTools Genotyping input channels
    channel_data = Channel.fromPath(params.channel_data_file, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')

    // BCFTools Genotyping input channels
    input_file_ids_bams_indexes = channel_data.map { row -> tuple(row.file_id, row.bam_file, row.bam_index_file) }
    file_id_reference_files_ch = channel_data.map { row -> tuple(row.file_id, row.sample_key, row.reference_file, row.snp_list) }

    // Run BCFTools Genotyping workflow
    GENOTYPING_BCFTOOLS(input_file_ids_bams_indexes, file_id_reference_files_ch)
}