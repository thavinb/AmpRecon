#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { DESIGNATE_PANEL_RESOURCES } from './designate_panel_resources.nf'
include { PULL_FROM_IRODS } from './pipeline-subworkflows/pull_from_irods.nf'
/*
include { validate_general_params } from '../main.nf'

// check paramater functions definition ------------------------------------
def validate_irods_exclusive_params(){
*/
    /*
    This functions counts the number of errors on input parameters exclusively used on IRODS subworkflow
    
    checks:
     - if irods manifest was provided
     - if irods manifest provided exists
    
    False for any of those conditions counts as an error.
    
    Returns
    ---
    <int> the nunber of errors found
    */
 /*   
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
    validate_general_params()
    // count errors and kill nextflow if any had been found
    erros += validate_irods_exclusive_params(errors)

    // check IRODS params
    if (params.execution_mode == "irods"){
        erros += validate_irods_exclusive_params(errors)
    }

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}
*/
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
        irods_ch.map{it -> tuple(it[0], it[1])}.set{new_sample_tag_panel_ch} // tuple (new_sample_id, panel_name)
        DESIGNATE_PANEL_RESOURCES(new_sample_tag_panel_ch, reference_ch)
        sample_tag_reference_files_ch = DESIGNATE_PANEL_RESOURCES.out.sample_tag_reference_files_ch
        // tuple (new_sample_id, panel_name, path/to/reference/genome, snp_list)

        // run step1.2b - pull from iRODS
        PULL_FROM_IRODS(irods_ch.map{it -> tuple(it[0], it[2])}) // tuple (new_sample_id, irods_path)
        bam_files_ch = PULL_FROM_IRODS.out.bam_files_ch

    emit:
        bam_files_ch
        sample_tag_reference_files_ch // tuple (new_sample_id, panel_name, path/to/reference/genome, snp_list)
	    file_id_to_sample_id_ch // tuple (file_id, sample_id)
}

