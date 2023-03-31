#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import irods processes
include { irods_retrieve } from '../modules/irods_retrieve.nf'
include { scramble_cram_to_bam } from '../modules/scramble.nf'

include { validate_general_params } from '../main.nf'

// check paramater functions definition ------------------------------------
def count_irods_to_reads_params_errors(){

    /*
    This functions counts the number of errors on input parameters exclusively used on IRODS subworkflow
    
    checks:
     - if irods manifest was provided
     - if irods manifest provided exists
    
    False for any of those conditions counts as an error.
    
    Returns
    ---
    <int> the number of errors found
    */
    def err = 0
    if (params.irods_manifest == null){
        log.error("An irods_manifest parameter must be provided for execution mode '${params.execution_mode}'.")
        err += 1
    }

    if (params.irods_manifest){
        irods_manifest = file(params.irods_manifest)
        if (!irods_manifest.exists()){
            log.error("The irods manifest file specified (${params.irods_manifest}) does not exist.")
            err += 1
        }
        else {
            validate_irods_mnf(params.irods_manifest, params.panels_settings)
        }
    }

    return err
}

def validate_parameters(){
    def errors = 0
    // import general params check ones
    errors += validate_general_params()
    // count errors and kill nextflow if any had been found
    errors += validate_irods_exclusive_params()

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

// -------------------------------------------------------------------------

workflow SANGER_IRODS_TO_READS {
    take:
        irods_manifest // irods manifest file
        reference_ch // tuple (fasta, panel_name, snp_list)
    main:
        // load manifest content
        irods_ch =  Channel.fromPath(irods_manifest, checkIfExists: true)
                      | splitCsv(header: true, sep: '\t')
                      //| map { row -> tuple(row.id_run, row.primer_panel, row.WG_lane) }
                      | map { row ->
                              WG_lane = "${row.irods_path}".split('/')[-1].split('\\.')[0]
                              tuple(row.sample_id, row.primer_panel, WG_lane, row.irods_path) 
                            }
                      | map { it -> tuple("${it[2]}_${it[0]}_${it[1]}", it[1], it[3], it[0])} // tuple (WG_lane_sample_id_panel_name, panel_name, irods_path, sample_id)

        // link file_id to sample_id
        irods_ch.map{it -> tuple(it[0], it[3])}.set{file_id_to_sample_id_ch}

        // assign each sample tag the appropriate set of reference files
        irods_ch.map{it -> tuple(it[0], it[1])}.set{new_file_id_panel_ch} // tuple (file_id, panel_name)

        new_file_id_panel_ch
          |  combine(reference_ch,  by: 1) // tuple (panel_name, file_id, fasta,snp_list)
          |  map{it -> tuple(it[1], it[0], it[2], it[3])}
          |  set{file_id_reference_files_ch}
        // tuple (file_id, panel_name, path/to/reference/genome, snp_list)

        // Retrieve CRAM files from iRODS
        irods_retrieve(irods_ch)

        // Convert iRODS CRAM files to BAM format
        scramble_cram_to_bam(irods_retrieve.out)
        bam_files_ch = scramble_cram_to_bam.out

    emit:
        bam_files_ch // tuple(file_id, bam_file)
        file_id_reference_files_ch // tuple (new_sample_id, panel_name, path/to/reference/genome, snp_list)
	    file_id_to_sample_id_ch // tuple (file_id, sample_id)
}

