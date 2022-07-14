#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../modules/bam_reset.nf'
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'
include { align_bam } from '../modules/align_bam.nf'

workflow redo_alignment {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

  take:
    sample_tag
    bam_file
    run_id
    reference_fasta_new
    reference_idx_fls
  main:
    // Unmap the bam files (ubam)
    bam_reset(sample_tag, bam_file)

    // convert ubams to fastqs
    bam_to_fastq(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)

    // do new alignment
    align_bam(
      bam_to_fastq.out.sample_tag,
      bam_to_fastq.out.fastq,
      reference_fasta_new,
      reference_idx_fls,
      run_id
      )
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
