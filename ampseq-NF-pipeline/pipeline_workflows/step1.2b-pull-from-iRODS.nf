#!/usr/bin/env nextflow
// enable dsl2
nextflow.enable.dsl = 2

// import irods processes
include { irods_manifest_parser } from '../modules/irods_manifest_parser.nf'
include { irods_retrieve } from '../modules/irods_retrieve.nf'
include { scramble_cram_to_bam } from '../modules/scramble.nf'

process writeOutputManifest {

  //publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_tag), path(bam_file), val(run_id)
    // TODO create a python box container

  output:
    //tuple val(run_id), path("${params.results_dir}/${run_id}/${run_id}_out1.2_mnf.csv")
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

workflow pull_from_iRODS {
  take:
    irods_ch
    ref_fasta
    ref_fasta_index_fl

  main:
    // process manifest
    // Parse iRODS manifest file
    irods_manifest_parser(irods_ch)

    // Retrieve CRAM files from iRODS
    irods_retrieve(irods_manifest_parser.out)

    // Convert iRODS CRAM files to BAM format
    scramble_cram_to_bam(irods_retrieve.out,
                         ref_fasta,
                         ref_fasta_index_fl)

    // Concatenate in-country BAM channel with iRODS BAM channel
    bam_files_ch = scramble_cram_to_bam.out
    //bam_ch.concat(scramble_cram_to_bam.out).set{ bam_files_ch }

    // Write manifest 1.2b_out.csv
    writeOutputManifest(bam_files_ch)
    //writeOutputManifest.out.view()
    bam_files_ch.view()
    // --------------------------------------------------------------------

  emit:
  bam_files_ch
}

/*
// -------------------------- DOCUMENTATION -----------------------------------
*/
