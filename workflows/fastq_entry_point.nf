#!/usr/bin/env nextflow

/*
    | FASTQ_ENTRY_POINT |-----------------------------------------
    
    This workflow takes a manifest file containing the paths to the 
    FASTQ files. 

    ------------------------------------------------------------------
*/

// enable dsl2
nextflow.enable.dsl = 2

// import processes
include { parse_panel_settings } from '../modules/parse_panels_settings.nf'
include { validate_fastq_mnf } from '../modules/validate_fastq_mnf.nf'

// -------------------------------------------------------------------------

workflow FASTQ_ENTRY_POINT {
    take:
        fastq_manifest // fastq manifest file
        reference_ch // tuple (reference_fasta_file, panel_name, snp_list)
    main:        
        index=0
        // get the sample data from the manifest and create file_id (combination of sample_id and panel name)
        fastq_ch = fastq_manifest
                    | splitCsv(header: true,sep: '\t')
                    | map { row ->
                            //added index to the file_id 
                            index = index+1
                            tuple(row.sample_id, row.primer_panel, row.fastq_path, "${row.sample_id}_${row.primer_panel}_${index}")
    }
        //link everything from fastq_manifest with reference_ch
        fastq_ch
          |  combine(reference_ch,  by: 1) // tuple (panel_name, sample_id, fastq_file, file_id, reference_fasta_file, snp_list)
          |  set{fastq_ch_comb_reference_ch}
        
        //create fastq_files_ch: file_id, fastq_file, reference_fasta_file
        fastq_ch_comb_reference_ch
          |  map{it -> tuple(it[3], it[2], it[4])} // tuple (file_id, fastq_file, reference_fasta_file)
          |  set{fastq_files_ch}

        //create file_id_reference_files_ch: file_id, panel_name, reference_fasta_file, snp_list
        fastq_ch_comb_reference_ch
          |  map{it -> tuple(it[3], it[0], it[4], it[5])} //tuple (file_id, panel_name, reference_fasta_file, snp_list)
          |  set{file_id_reference_files_ch}

        // link 'file_id' to sample_id
        fastq_ch.map{it -> tuple(it[3],it[0])}.set{file_id_to_sample_id_ch}    // tuple(file_id, sample_id)

    emit:
        fastq_files_ch // tuple (file_id, fastq_file, reference_fasta_file)
        file_id_reference_files_ch // tuple (file_id, panel_name, reference_fasta_file, snp_list)
        file_id_to_sample_id_ch // tuple (file_id, sample_id)
}

def fastq_parameter_check(){

    /*
    This functions counts the number of errors on input parameters exclusively used on FASTQ entry point subworkflow
    
    checks:
     - if fastq manifest was provided
     - if fastq manifest provided exists
    
    False for any of those conditions counts as an error.
    
    Returns
    ---
    <int> the number of errors found
    */

    def error = 0
    if (params.fastq_manifest == null){
        log.error("An fastq_manifest parameter must be provided for execution mode '${params.execution_mode}'.")
        error += 1
    }

    if (params.fastq_manifest){
        fastq_manifest = file(params.fastq_manifest)
        if (!fastq_manifest.exists()){
            log.error("The fastq manifest file specified (${params.fastq_manifest}) does not exist.")
            error += 1
        }
        else {
            validate_fastq_mnf(params.fastq_manifest, params.panels_settings)
        }
    }

    if (error > 0) {
        log.error("Parameter errors were found, the pipeline will not run")
        exit 1
    }
}


workflow {
    
    ref_and_annt_ch = parse_panel_settings(params.panels_settings)
    reference_ch = ref_and_annt_ch[0] // tuple(reference_file, panel_name, snp_list)
    
    fastq_manifest = Channel.fromPath(params.fastq_manifest)
    
    fastq_parameter_check()

    // Run FASTQ_ENTRY_POINT workflow
    FASTQ_ENTRY_POINT(fastq_manifest, reference_ch)
}
