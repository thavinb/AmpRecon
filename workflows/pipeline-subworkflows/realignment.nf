#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../../modules/bam_reset.nf'
include { bam_to_fastq } from '../../modules/bam_to_fastq.nf'
include { align_bam } from '../../modules/align_bam.nf'
include { scramble_sam_to_bam } from '../../modules/scramble.nf'
include { samtools_sort } from '../../modules/samtools.nf'
include { samtools_index } from '../../modules/samtools.nf'
<<<<<<< HEAD
include { files_and_panels_to_csv } from '../../modules/read_counts_per_region.nf'
include { read_count_per_region } from '../../modules/read_counts_per_region.nf'    
=======
include { bam_ref_ch_to_csv } from '../../modules/read_counts_per_region.nf'
include { read_count_per_region } from '../../modules/read_counts_per_region.nf'    
include { files_and_panels_to_csv } from '../../modules/read_counts_per_region.nf'
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
include { upload_pipeline_output_to_s3 } from '../../modules/upload_pipeline_output_to_s3.nf'

workflow REALIGNMENT {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

  take:
    sample_tag
    bam_file
    run_id
    sample_tag_reference_files_ch // tuple (sample_id, fasta_file, [fasta_indx_files], panel_name)
<<<<<<< HEAD
    pannel_anotations_files // tuple (pannel_name, anotation_file)


=======
    annotations_ch // tuple (pannel_name, anotation_file)

>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
  main:
    // Unmap the bam files (ubam)
    bam_reset(sample_tag, bam_file)
    
    // convert ubams to fastqs
<<<<<<< HEAD
    bam_to_fastq(bam_reset.out)
=======
    bam_to_fastq(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
    // tuple (sample_tag, fastq)

    // prepare channels to be used on join for input for other processes

    bam_to_fastq.out // tuple (sample_id, fastq_file)
          | join(sample_tag_reference_files_ch) //tuple (sample_id, fastq, fasta_file, fasta_idx_files, panel_name)
          | set{align_bam_In_ch}

     // do new alignment
    align_bam(align_bam_In_ch)

    // convert sam to bam
    scramble_sam_to_bam(align_bam.out.sample_tag, align_bam.out.sam_file)
    
    // sort and index bam
    samtools_sort(scramble_sam_to_bam.out)
    samtools_index(samtools_sort.out)

    // DO READCOUNTS 

    // make CSV file of bam file names with associated amplicon panel
    bam_file_names = samtools_index.out.map{it -> it[1].baseName}.collect() // amplicon panel now part of BAM name
    files_and_panels_to_csv(bam_file_names)
    bams_and_indices = samtools_index.out.map{it -> tuple(it[1], it[2])}.collect()

    read_count_per_region(
        run_id,
        files_and_panels_to_csv.out,
        bams_and_indices,
<<<<<<< HEAD
        pannel_anotations_files
=======
        annotations_ch,
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
    )
    // upload read counts and BAM files to S3 bucket
    output_bams_ch = samtools_index.out.map{it -> it[1]}

    read_count_per_region.out.qc_csv_file.concat(output_bams_ch).set{output_to_s3}
    upload_pipeline_output_to_s3(output_to_s3)

  //emit:

}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Collated BAM files are reset to their prealigned state by removing all
     @SQ from header, all reads marked as unmapped, dropping non-primary
     alignments and sorting order to set to unknown.
[2] BAM files are converted to FASTQ format, before being aligned to a reference genome.
[3] align to a given reference
*/


