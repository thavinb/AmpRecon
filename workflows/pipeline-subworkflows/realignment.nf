#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../../modules/bam_reset.nf'
include { bam_to_fastq } from '../../modules/bam_to_fastq.nf'
include { align_bam } from '../../modules/align_bam.nf'
include { scramble_sam_to_bam } from '../../modules/scramble.nf'
include { add_read_group } from '../../modules/add_read_group.nf'
include { samtools_sort } from '../../modules/samtools.nf'
include { samtools_index } from '../../modules/samtools.nf'
include { files_and_panels_to_csv } from '../../modules/read_counts_per_region.nf'
include { read_count_per_region } from '../../modules/read_counts_per_region.nf'
include { upload_pipeline_output_to_s3 } from '../../modules/upload_pipeline_output_to_s3.nf'

workflow REALIGNMENT {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

  take:
    sample_tag
    bam_file
    sample_tag_reference_files_ch // tuple (sample_id, fasta_file, [fasta_indx_files], panel_name)
    annotations_ch // tuple (pannel_name, anotation_file)

  main:
    // Unmap the bam files (ubam)
    bam_reset(sample_tag, bam_file)
    
    // convert ubams to fastqs
    bam_to_fastq(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)

    // prepare channels to be used on join for input for other processes
    bam_to_fastq.out // tuple (sample_id, fastq_file)
          | join(sample_tag_reference_files_ch) //tuple (sample_id, fastq, fasta_file, fasta_idx_files, panel_name)
          | set{align_bam_In_ch}

     // do new alignment
    align_bam(align_bam_In_ch)

    // convert sam to bam
    scramble_sam_to_bam(align_bam.out.sample_tag, align_bam.out.sam_file)

    // add correct read group from reset bams into aligned bams
    scramble_sam_to_bam.out.join(bam_reset.out.bam_reset_tuple_ch).set{add_read_group_input_ch}
    add_read_group(add_read_group_input_ch)
    
    // sort and index bam
    samtools_sort(add_read_group.out)
    samtools_index(samtools_sort.out)

    // DO READCOUNTS 

    // make CSV file of bam file names with associated amplicon panel
    bam_file_names = samtools_index.out.map{it -> it[1].baseName}.collect() // amplicon panel now part of BAM name
    files_and_panels_to_csv(bam_file_names)
    bams_and_indices = samtools_index.out.map{it -> tuple(it[1], it[2])}.collect()

    read_count_per_region(
        files_and_panels_to_csv.out,
        bams_and_indices,
        annotations_ch,
    )

    // upload read counts and BAM files / indices to S3 bucket
    if (params.upload_to_s3){
      output_bams_and_indices_ch = samtools_index.out.map{it -> tuple(it[1], it[2])}.flatten()
      read_count_per_region.out.qc_csv_file.concat(output_bams_and_indices_ch).set{output_to_s3}
      upload_pipeline_output_to_s3(output_to_s3)
    }

  emit:
    samtools_index.out
}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Collated BAM files are reset to their prealigned state by removing all
     @SQ from header, all reads marked as unmapped, dropping non-primary
     alignments and sorting order to set to unknown.
[2] BAM files are converted to FASTQ format, before being aligned to a reference genome.
[3] align to a given reference
*/


