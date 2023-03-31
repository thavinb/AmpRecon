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

process write_aligned_bam_mnf {

  //publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true
  label 'pythonBox'
  input:
    tuple val(file_id), path(bam_file), path(bam_idx), val(panel_name)

// The $/ ... /$ is necessary to avoid nextflow to read "\n" incorrectly
$/
#!/usr/bin/python3
from pathlib import Path

# setup inputs
bam_fl = "${bam_file}"
bam_idx = "${bam_idx}"
file_id = "${file_id}"
publishDir = f"${params.results_dir}/"
bam_dir=f"${params.results_dir}/"
panel_name = f"${panel_name}"
# if manifest already exists, just append new lines
path_to_mnf = f"{publishDir}/aligned_bams_mnf.csv"
if Path(path_to_mnf).is_file():
    out_mnf = open(f"{path_to_mnf}", "a")

# if manifest does not exist, create file and write header
else:
    out_mnf = open(f"{path_to_mnf}", "w")
    out_mnf.write("file_id,panel_name,bam_file,bam_idx\n")

# write manifest line for the bam file
out_mnf.write(f"{file_id},{panel_name},{bam_dir}{bam_fl},{bam_dir}{bam_idx}\n")
out_mnf.close()
/$
}

workflow ALIGNMENT {
  //remove alignment from bam - this process proceeds directly after the end of 1.2x

  take:
    file_id
    bam_file
    file_id_reference_files_ch // tuple (file_id, fasta_file, panel_name)

  main:
    // Unmap the bam files (ubam)
    bam_reset(file_id, bam_file)
    
    // convert ubams to fastqs
    bam_to_fastq(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)

    // prepare channels to be used on join for input for other processes
    bam_to_fastq.out // tuple (file_id, fastq_file)
          | join(file_id_reference_files_ch) //tuple (file_id, fastq, fasta_file, panel_name)
          | set{ bwa_ch }

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
  // write aligned bam manifest

  samtools_index.out // tuple (file_id, input_bam, bam_bai)
    | join(file_id_reference_files_ch) //tuple (file_id, input_bam, bam_bai, fasta_file, panel_name)
    | map {it -> tuple(it[0], it[1], it[2], it[4])}
    | set {manifest_ch} // tuple (file_id, bam, bam_index, panel_name)

  write_aligned_bam_mnf(manifest_ch)

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
