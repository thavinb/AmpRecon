#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { collate_alignments } from './modules/collate_alignments.nf'
include { bam_reset } from './modules/bam_reset.nf'
include { clip_adapters } from './modules/clip_adapters.nf'
include { bam_to_fastq } from './modules/bam_to_fastq.nf'
include { align_bam } from './modules/align_bam.nf'
include { scramble_sam_to_bam } from './modules/scramble.nf'
include { mapping_reheader } from './modules/mapping_reheader.nf'
include { bam_split } from './modules/bam_split.nf'
include { bam_merge } from './modules/bam_merge.nf'
include { alignment_filter } from './modules/alignment_filter.nf'
include { sort_bam } from './modules/sort_bam.nf'

def load_intermediate_ch(csv_ch){
  // if not set as parameter, assumes is a channel containing a path for the csv
  intermediateCSV_ch = csv_ch |
                splitCsv(header:true) |
                multiMap {row -> run_id:row.run_id
                                 cram_fl:row.cram_fl
                                 sample_tag:row.sample_tag}
  return intermediateCSV_ch
}

process writeOutputManifest {

  //publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_tag), path(bam_file)
    val(run_id)
    // TODO create a python box container

  //output:
  //  tuple val(run_id), path("${run_id}_out1.2_mnf.csv")
// The $/ ... /$ is necessary to avoid nextflow to read "\n" incorrectly
$/
#!/usr/bin/python3
from pathlib import Path

# setup inputs
run_id = "${run_id}"
bam_fl = "${bam_file}"
sample_tag = "${sample_tag}"
publishDir = f"${params.results_dir}/{run_id}/"
bam_dir=f"${params.results_dir}{run_id}/"

# if manifest already exists, just append new lines
path_to_mnf = f"{publishDir}/{run_id}_out1.2_mnf.csv"
if Path(path_to_mnf).is_file():
    out_mnf = open(f"{path_to_mnf}", "a")

# if manifest does not exist, create file and write header
else:
    out_mnf = open(f"{path_to_mnf}", "w")
    out_mnf.write("run_id,bam_fl,sample_tag\n")

# write manifest line for the bam file
out_mnf.write(f"{run_id},{bam_dir}{bam_fl},{sample_tag}\n")
out_mnf.close()
/$
}

workflow cram_to_bam {
    take:
        // manifest from step 1.1
        intermediate_csv
        sample_tag_reference_files_ch

    main:
        // Process manifest
        intCSV_ch = load_intermediate_ch(intermediate_csv)

        // Collate cram files by name
        collate_alignments(intCSV_ch.run_id,
                          intCSV_ch.cram_fl,
                          intCSV_ch.sample_tag)

        // Transform BAM file to pre-aligned state
        bam_reset(collate_alignments.out.sample_tag,
                  collate_alignments.out.collated_bam)

        // Remove adapters
        clip_adapters(
                  bam_reset.out.sample_tag,
                  bam_reset.out.reset_bam)

        // Convert BAM to FASTQ
        bam_to_fastq(clip_adapters.out.tuple)

        bam_to_fastq.out.join(sample_tag_reference_files_ch).combine(intCSV_ch.run_id).unique().set{align_bam_In_ch}
        //align_bam_In_ch = [sample_tag, fastq, ref_fasta, ref_bwa_idx_fls, run_id]
        //align_bam_In_ch.view()
        align_bam(align_bam_In_ch)

        // Convert SAM to BAM
        //bambi_select(align_bam.out.sample_tag, align_bam.out.sam_file)
        scramble_sam_to_bam(align_bam.out.sample_tag,
                            align_bam.out.sam_file,
        )

        // Merges the current headers with the old ones.
        // Keeping @SQ.*\tSN:([^\t]+) matching lines from the new header.
        reheader_in_ch = scramble_sam_to_bam.out
                            | join(clip_adapters.out.tuple)
                            | join(sample_tag_reference_files_ch)

        mapping_reheader(reheader_in_ch)

        // Split BAM rank pairs to single ranks per read
        bam_split(mapping_reheader.out)

        // Merge BAM files with same reads
        bam_merge_In_ch = bam_split.out.join(clip_adapters.out.tuple)

        bam_merge(bam_merge_In_ch)

        // Split alignments into different files
        alignment_filter(bam_merge.out.sample_tag,
                         bam_merge.out.merged_bam)

        // BAM sort by coordinate
        sort_bam(intCSV_ch.run_id,
                 alignment_filter.out.sample_tag,
                 alignment_filter.out.selected_bam)
        bam_ch = sort_bam.out

        // --------------------------------------------------------------------
        // write manifest out
        writeOutputManifest(bam_ch, intCSV_ch.run_id)

    emit:
        bam_ch
}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Read pairs within each CRAM file are first collated and output in BAM format.
[2] Collated BAM files are reset to their prealigned state by removing all @SQ from header, all reads marked as unmapped,
dropping non-primary alignments and sorting order to set to unknown.
[4] Previously identified adapter sequences in each BAM are then clipped off and moved into auxiliary fields.
[5] BAM files are converted to FASTQ format, before being aligned to a reference genome.
[6] The mapped SAM files are scrambled / converted into BAM files.
[7] Specific headers from the adapters clipped BAM file are copied to the scrambled, realigned BAM files.
Duplicate IDs for @RG and @PG records in the header are made unique with the addition of a suffix, read records are also updated.
[8] Rank pairs produced by the collation are converted into single ranks per read.
[9] BAM files with the same reads are merged. Ranks are also stripped from the alignments, clipped sequences are reinserted and quality string parts added.
[10] Reads from all of the BAM files are split into separate files depending on their alignments.
[11] Alignments in the BAM files are sorted by coordinate. The sorted BAM files are emitted ready for further processing.

*/

