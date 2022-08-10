#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from './modules/bam_reset.nf'
include { bam_to_fastq } from './modules/bam_to_fastq.nf'
include { align_bam } from './modules/align_bam.nf'
include { scramble_sam_to_bam } from './modules/scramble.nf'
include { sort_and_index } from './modules/read_count_per_region.nf'
include { read_count_per_region_qc } from './modules/read_count_per_region.nf'
include { get_sample_ref } from '../step1.2a-cram-to-bam/modules/get_sample_ref.nf'

def load_intermediate_ch(csv_ch){
  // if not set as parameter, assumes is a channel containing a path for the csv
  intermediateCSV_ch = csv_ch |
                splitCsv(header:true) |
                multiMap {row -> run_id:row.run_id
                                 cram_fl:row.cram_fl
                                 sample_tag:row.sample_tag}
  return intermediateCSV_ch
}

workflow redo_alignment {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

  take:
    intermediate_csv

  main:
    // Process manifest
    intCSV_ch = load_intermediate_ch(intermediate_csv)
    //--SAMPLE TO REF RELATIONSHIP ----------------------------------------
    // get sample references from [run_id]_manifest.csv
    // curent solution assumes [run_id]_manifest.csv is at the output folder
    // TODO: in the future, changed it to be specified
    get_sample_ref(
        intCSV_ch.run_id,
        intCSV_ch.sample_tag,
        intCSV_ch.bam_fl
    )
    sample_ref_ch = get_sample_ref.out

    // prepare channels to be used on join for input for other processes
    sample_ref_ch.map {it -> tuple(it[2],it[3],it[4], it[0])}.set{sample2ref_tuple_ch}

    // Unmap the bam files (ubam)
    bam_reset(intCSV_ch.sample_tag, intCSV_ch.bam_fl)

    // convert ubams to fastqs
    bam_to_fastq(bam_reset.out.sample_tag,
        bam_reset.out.reset_bam)

    bam_to_fastq.out.join(sample2ref_tuple_ch).set{align_bam_In_ch}

    // do new alignment
    align_bam(align_bam_In_ch)

    // convert sam to bam_dir
    scramble_sam_to_bam(align_bam.out.sample_tag, align_bam.out.sam_file)

    // sort and index bam
    sort_and_index(scramble_sam_to_bam.out)
    sort_and_index.out.bam_dir.unique().collect().set{bam_dir_ch} // Needed to ensure correct execution order

    // Get read counts
    qc_run_ids_ch = Channel.from("GRC1", "GRC2", "Spec")
    qc_run_cnf_files_ch = Channel.from(file(params.grc1_qc_file), file(params.grc2_qc_file), file(params.spec_qc_file))
    read_count_per_region_qc(
        run_id,
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
