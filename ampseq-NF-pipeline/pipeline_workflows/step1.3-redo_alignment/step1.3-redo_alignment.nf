#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from './modules/bam_reset.nf'
include { bam_to_fastq } from './modules/bam_to_fastq.nf'
include { align_bam } from './modules/align_bam.nf'
include { scramble_sam_to_bam } from './modules/scramble.nf'
include { samtools_sort } from './modules/samtools.nf'
include { samtools_index } from './modules/samtools.nf'
include { bam_ref_ch_to_csv } from './modules/read_count_per_region.nf'
include { read_count_per_region } from './modules/read_count_per_region.nf'

workflow redo_alignment {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

take:
    sample_tag
    bam_file
    run_id
    sample_tag_reference_files_ch // tuple (sample_id, fasta_file, [fasta_indx_files])
    
  main:
    // Unmap the bam files (ubam)
    bam_reset(sample_tag, bam_file)
    
    // convert ubams to fastqs
    bam_to_fastq(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)
    // tuple (sample_tag, fastq)

    // prepare channels to be used on join for input for other processes
    
    bam_to_fastq.out // tuple (sample_id, fastq_file)
          | join(sample_tag_reference_files_ch) //tuple (sample_id, fastq, fasta_file, fasta_idx_files)
          | set{align_bam_In_ch}
 
     // do new alignment
    align_bam(align_bam_In_ch)

    // convert sam to bam_dir
    scramble_sam_to_bam(align_bam.out.sample_tag, align_bam.out.sam_file)
    
    // PREPARATION FOR READ COUNTS 
    // sort and index bam
    samtools_sort(scramble_sam_to_bam.out)
    samtools_index(samtools_sort.out.bam)
    
    // GAMBIARRA ALERT --------------------------------------------
    // Needed to ensure correct execution order
    samtools_index.out.bam_dir.unique().collect().set{bam_dir_ch} 
    // ------------------------------------------------------------
     
    // join references to indexed bam channel
    samtools_index.out.files // tuple (sample_tag, input_bam, input_bai)
                | join(sample_tag_reference_files_ch) // (sample_tag, input_bam, input_bai, fasta_file, fasta_idx_file)
                | map {it -> tuple(it[1].simpleName, it[3])} // (input_bam_name,fasta_file)
                | set{bam_ref_ch}

    // output channel to csv - used for making read count plex files
    bam_ref_ch_to_csv(bam_ref_ch)

    // get read counts
    qc_run_ids_ch = Channel.from("GRC1", "GRC2", "Spec")
    qc_run_cnf_files_ch = Channel.from(
                          file(params.grc1_qc_file),
                          file(params.grc2_qc_file),
                          file(params.spec_qc_file))

    read_count_per_region(
        run_id,
        "${launchDir}/bam_ref_ch.csv",
        bam_dir_ch,
        qc_run_ids_ch,
        qc_run_cnf_files_ch
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
