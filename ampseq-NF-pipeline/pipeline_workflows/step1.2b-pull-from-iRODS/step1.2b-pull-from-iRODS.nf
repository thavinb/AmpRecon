#!/usr/bin/env nextflow
// enable dsl2
nextflow.enable.dsl = 2

// import irods processes
include { irods_manifest_parser } from './modules/irods_manifest_parser.nf'
include { irods_retrieve } from './modules/irods_retrieve.nf'
include { scramble_cram_to_bam } from './modules/scramble.nf'

process writeOutputManifest {

  input:
    tuple val(sample_tag), path(bam_file), val(run_id)

  output:
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
    irods_ch // tuple(id_run, WG_lane)
    sample_id_ref_ch // tuple(WG_lane, primer_panel, fasta_files) * needed to associate new_ids to pannels
  main:
    // get names and paths 
    irods_manifest_parser(irods_ch)
    //remove WG_lane (not needed on other processes), but must be associated with "new" sample id
    //    tuple new_sample_id, ${iRODS_file_path}, id_run, WG_lane

    irods_manifest_parser.out.map{it -> tuple( it[0], it[1], it[2] ) }.set {irods_retrieve_In_ch}

    // get new_sample_ids pannels

    irods_manifest_parser.out //  tuple(new_sample_id, iRODS_file_path, id_run, WG_lane)
              | map {it -> tuple( it[3], it[0])}
              | set {newSample_WgLn_ch}  // tuple(WG_lane, new_sample_id)

    sample_id_ref_ch
              | map { it -> tuple( it[0], it[2])} // tuple(WG_lane, fasta_files)
              | join (newSample_WgLn_ch)  // tuple(WG_lane, fasta_files, new_sample_id)
              | map { it -> tuple(it[2], it[1])} // tuple( new_sample_id, fasta_files)
              | set { sample_to_ref_ch }
    
    // Retrieve CRAM files from iRODS
    irods_retrieve(irods_retrieve_In_ch)
    // Convert iRODS CRAM files to BAM format
    scramble_cram_to_bam(irods_retrieve.out)

    // Concatenate in-country BAM channel with iRODS BAM channel
    bam_files_ch = scramble_cram_to_bam.out
    //bam_ch.concat(scramble_cram_to_bam.out).set{ bam_files_ch }
    
 
    // Write manifest 1.2b_out.csv
    writeOutputManifest(bam_files_ch)

    // --------------------------------------------------------------------

  emit:
    bam_files_ch  // tuple(new_sample_id, bam_file)
    sample_to_ref_ch // tuple(new_sample_id, ref_files)
}

/*
// -------------------------- DOCUMENTATION -----------------------------------
*/
