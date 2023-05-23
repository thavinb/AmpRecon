#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bwa_alignment_and_post_processing } from '../modules/bwa_alignment_and_post_processing.nf'
include { upload_pipeline_output_to_s3 } from '../modules/upload_pipeline_output_to_s3.nf'

workflow ALIGNMENT {
  take:
    fastq_ch // tuple (file_id, fastq_file, reference_fasta_file)

  main:
    // do alignment, sort by coordinate and index
    bwa_alignment_and_post_processing(fastq_ch)

    // upload BAM files and index files to S3 bucket
    if (params.upload_to_s3){
      output_to_s3 = bwa_alignment_and_post_processing.out.map{it -> tuple(it[1], it[2])}.flatten()
      upload_pipeline_output_to_s3(output_to_s3)
    }

  emit:
    bwa_alignment_and_post_processing.out // tuple (file_id, bam_file, bai_file)
}

