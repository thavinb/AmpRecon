#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { collate_alignments } from '../modules/collate_alignments.nf'
include { bam_reset } from '../modules/bam_reset.nf'
include { clip_adapters } from '../modules/clip_adapters.nf'
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'
include { align_bam } from '../modules/align_bam.nf'
include { bambi_select  } from '../modules/scramble_sam_to_bam.nf'
include { mapping_reheader } from '../modules/mapping_reheader.nf'
include { bam_split } from '../modules/bam_split.nf'
include { bam_merge } from '../modules/bam_merge.nf'
include { alignment_filter } from '../modules/alignment_filter.nf'
//include { sort_bam } from '../modules/sort_bam.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'

def load_manifest_ch(csv_ch){
  //if csv file is provided as parameter, use it by default and ignore input
  if (!(params.manifest_step1_2 == '')){
      // TODO : add check if file exist
      manifest_fl = params.manifest_step1_2
      csv_ch = Channel.fromPath(manifest_fl)
      }
  // if not set as parameter, assumes is a channel containing a path for the csv
  manifest_ch = csv_ch |
                splitCsv(header:true) |
                multiMap {row -> run_id:row.run_id
                                 cram_fl:row.cram_fl
                                 sample_tag:row.sample_tag}
                //map {row-> tuple(val(row.run_id), path(row.cram_fl), val(row.tag))}

  return manifest_ch
}

workflow cram_to_bam {
    take:
        // manifest from step 1.1
        manifest_fl
        // reference index files
        ref_bwa_index_fls
        ref_fasta_index_fl
        ref_dict_fl
    main:
        // Process manifest
        mnf_ch = load_manifest_ch(manifest_fl)

        // Collate cram files by name
        collate_alignments(mnf_ch.run_id, mnf_ch.cram_fl, mnf_ch.sample_tag)
        //collate_alignments.out.view()

        // Transform BAM file to pre-aligned state

        bam_reset(collate_alignments.out.sample_tag,
                  collate_alignments.out.collated_bam)
        bamReset_Out_ch = bam_reset.out
        // Remove adapters
        clip_adapters(bam_reset.out.sample_tag, bam_reset.out.prealigned_bam)

        // Convert BAM to FASTQ
        bam_to_fastq(clip_adapters.out.sample_tag,
                     clip_adapters.out.clipped_bam)

        align_bam(bam_to_fastq.out.sample_tag,
                  bam_to_fastq.out.fastq,
                  params.reference_fasta,
                  ref_bwa_index_fls)

        // Map FASTQ to reference

        // SAM to BAM
        // scramble sam to bam (?)
        bambi_select(align_bam.out.sample_tag, align_bam.out.sam_file)

        // Merges the current headers with the old ones. Keeping @SQ.*\tSN:([^\t]+) matching lines from the new header.
        /*
        scramble_sam_to_bam.out
        .join(scramble_sam_to_bam.out)
        .join(all_manifest_data)
        .multiMap {
            tag: it[0]
            clipped: it[2]
            bam: it[3]
            reference_filename: it[4]["reference_fasta"].getName()
        }.set { mapping_reheader_input }
        mapping_reheader(mapping_reheader_input.tag, mapping_reheader_input.clipped, mapping_reheader_input.bam, mapping_reheader_input.reference_file, reference_index_files)

        // Split BAM rank pairs to single ranks per read
        mapping_reheader.out
        .join(bam_to_fastq.out)
        .multiMap {
            tag: it[0]
            reheadered_bam: it[1]
            clipped: it[3]
        }.set { bam_split_input }

        bam_split(bam_split_input.tag, bam_split_input.reheadered_bam, bam_split_input.clipped)

        // Merge BAM files with same reads
        bam_merge(bam_split.out, bam_split_input.clipped)

        // Split alignments into different files
        bam_merge.out
        .multiMap {
            tag: it[0]
            bam: it[1]
        }.set { alignment_filter_input }
        alignment_filter(alignment_filter_input.tag, alignment_filter_input.bam)

        // BAM sort by coordinate
        alignment_filter.out
        .multiMap {
            tag: it[0]
            alignment_file: it[1]
        }.set { bam_sort_input }
        sort_bam(bam_sort_input.tag, bam_sort_input.alignment_file)
        bam_ch = sort_bam.out

    emit:
        bam_ch
*/
}
/*
// -------------------------- DOCUMENTATION ----------------------------------- 
[1] Read pairs within each CRAM file are first collated and output in BAM format.
[2] Collated BAM files are reset to their prealigned state by removing all @SQ from header, all reads marked as unmapped, 
dropping non-primary alignments and sorting order to set to unknown. 
[4] Previously identified adapter sequences in each BAM are then clipped off and moved into auxiliary fields.
[5] BAM files are converted to FASTQ format, before being aligned to a reference genome.
[6] Reads in the mapped SAM files are split into separate BAM files depending on their alignments.
[7] Specific headers from the adapters clipped BAM file are copied to the scrambled, realigned BAM. Duplicate IDs for @RG and @PG records in the header 
are made unique with the addition of a suffix, read records are also updated.
[8] Rank pairs produced by the collation are converted into single ranks per read. 
[9] BAM files with the same reads are merged. Ranks are also stripped from the alignments, clipped sequences are reinserted and quality string parts added.
[10] Reads from all of the BAM files are split into separate files depending on their alignments.
[11] Alignments in the BAM files are sorted by coordinate. The sorted BAM files are emitted ready for further processing.

*/
