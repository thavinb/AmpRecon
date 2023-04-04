#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../modules/bam_reset.nf'
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'
include { bwa_alignment } from '../modules/bwa_alignment.nf'
include { scramble_sam_to_bam } from '../modules/scramble.nf'
include { add_read_group } from '../modules/add_read_group.nf'
include { samtools_sort } from '../modules/samtools.nf'
include { samtools_index } from '../modules/samtools.nf'
include { upload_pipeline_output_to_s3 } from '../modules/upload_pipeline_output_to_s3.nf'

workflow ALIGNMENT {
  take:
    file_id
    bam_file
    file_id_reference_files_ch // tuple (file_id, fasta_file, panel_name)

  main:
     // do new alignment
    bwa_alignment(bwa_ch)

    // convert sam to bam
    scramble_sam_to_bam(bwa_alignment.out.sample_tag, bwa_alignment.out.sam_file)

    // add correct read group from reset bams into aligned bams
    scramble_sam_to_bam.out.join(bam_reset.out.bam_reset_tuple_ch).set{add_read_group_input_ch}
    add_read_group(add_read_group_input_ch)
    
    // sort and index bam
    samtools_sort(add_read_group.out)
    samtools_index(samtools_sort.out)

    // upload BAM files and index files to S3 bucket
    if (params.upload_to_s3){
      output_to_s3 = samtools_index.out.map{it -> tuple(it[1], it[2])}.flatten()
      upload_pipeline_output_to_s3(output_to_s3)
    }

  emit:
    samtools_index.out
}
