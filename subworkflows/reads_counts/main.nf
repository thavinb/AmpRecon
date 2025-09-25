#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

// enable dsl2
nextflow.enable.dsl = 2

include { files_and_panels_to_csv } from '../modules/read_counts_per_region/main'
include { read_count_per_region } from '../modules/read_counts_per_region.nf/main'

workflow READ_COUNTS {
  take:
    indexed_bams_ch
    file_id_reference_files_ch
    annotations_ch // tuple (panel_name, annotation_file)

  main:

    // make CSV file of bam file names with associated amplicon panel
    bam_file_names = indexed_bams_ch.map{it -> it[1].baseName}.collect() // amplicon panel now part of BAM name
    panel_names = indexed_bams_ch.join(file_id_reference_files_ch).map{it -> it[4]}.collect()
    files_and_panels_to_csv(bam_file_names, panel_names)
    bams_and_indices = indexed_bams_ch.map{it -> tuple(it[1], it[2])}.collect()

    // determine read counts per amplicon region
    read_count_per_region(
      files_and_panels_to_csv.out,
      bams_and_indices,
      annotations_ch,
    )

}

workflow {
    // File required for Read Counts input channels
    channel_data = Channel.fromPath(params.channel_data_file, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')

    // Read Counts input channels
    indexed_bams_ch = channel_data.map { row -> tuple(row.file_id, file(row.bam_file), file(row.bam_index_file)) }
    file_id_reference_files_ch = channel_data.map { row -> tuple(row.file_id, row.panel_name, row.panel_name, row.panel_name) }
    annotations_ch = channel_data.map { row -> tuple(row.panel_name, file(row.annotation_file)) }.unique()

    // Run Read Counts workflow
    READ_COUNTS(indexed_bams_ch, file_id_reference_files_ch, annotations_ch)
}
