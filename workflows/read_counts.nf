#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { files_and_panels_to_csv } from '../modules/read_counts_per_region.nf'
include { read_count_per_region } from '../modules/read_counts_per_region.nf'
include { upload_pipeline_output_to_s3 } from '../modules/upload_pipeline_output_to_s3.nf'

workflow READ_COUNTS {
  take:
    indexed_bams_ch
    annotations_ch // tuple (panel_name, anotation_file)

  main:

    // make CSV file of bam file names with associated amplicon panel
    bam_file_names = indexed_bams_ch.map{it -> it[1].baseName}.collect() // amplicon panel now part of BAM name
    panel_names = indexed_bams_ch.join(sample_tag_reference_files_ch).map{it -> it[4]}.collect()
    files_and_panels_to_csv(bam_file_names, panel_names)
    bams_and_indices = indexed_bams_ch.map{it -> tuple(it[1], it[2])}.collect()

    // determine read counts per amplicon region
    read_count_per_region(
        files_and_panels_to_csv.out,
        bams_and_indices,
        annotations_ch,
    )

    // upload read counts to S3 bucket
    if (params.upload_to_s3){
      read_count_per_region.out.set{output_to_s3}
      upload_pipeline_output_to_s3(output_to_s3)
    }
}