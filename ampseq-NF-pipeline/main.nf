#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { bcl_to_cram } from './pipeline_workflows/step1.1-bcl-to-cram.nf'
include { get_taglist_file } from './modules/manifest2tag.nf'
include { validate_parameters; load_manifest_ch } from './modules/inputHandling.nf'
include { make_samplesheet_manifest } from './modules/samplesheet_parser.nf'
include { validate_samplesheet_manifest } from './modules/samplesheet_manifest_validation.nf'
// logging info ----------------------------------------------------------------
// This part of the code is based on the FASTQC PIPELINE (https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf)

/*
* ANSI escape codes to color output messages, get date to use in results folder name
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

log.info """
        ===========================================
         AMPSEQ_0.0 (dev : prototype)
         Used parameters:
        -------------------------------------------
         --manifest            : ${params.manifest}
         --results_dir         : ${params.results_dir}
         Runtime data:
        -------------------------------------------
         Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
         Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
         Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
         Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
         ------------------------------------------
         """
         .stripIndent()

//Run container:          ${ANSI_GREEN}${workflow.container}${ANSI_RESET}


// logging info ----------------------------------------------------------------

/*
workflow set_up_params {
    main:
        Channel.fromPath("${params.bcl_dir_path}").set{bcls}
        Channel.value("${params.lane}").set{lane}
        Channel.value("${params.run_id}").set{run_id}
        Channel.value("${params.study_name}").set{study_name}
        Channel.fromPath("${params.manifest_file_path}").set{manifest}

	//manifest2tag
        get_tag_list_file(manifest, "library", "sample", study_name)
        barcodes_file = get_tag_list_file.out.taglist_file

    emit:
        bcls
        lane
        run_id
        study_name
        barcodes_file
}
*/

// Main entry-point workflow



workflow {
    // -- Pre Processing ------------------------------------------------------
    // validate input
    validate_parameters()

    // process manifest input
    manifest_ch = load_manifest_ch()

    // process samplesheets manifest (necessary to get barcodes) and validate it
    make_samplesheet_manifest(manifest_ch)
    validate_samplesheet_manifest(make_samplesheet_manifest.out)

    get_taglist_file_In_ch = manifest_ch.join(validate_samplesheet_manifest.out)

    get_taglist_file(get_taglist_file_In_ch)
    tag_files_ch = get_taglist_file.out

    step1_Input_ch = manifest_ch.join(tag_files_ch)

    // use generate barcode files (run_id.taglist) from samplesheet manifest
    //manifest2tag_In_ch = manifest_ch.join(ssht_manifest_ch)
    //get_taglist_file(ssht_manifest_ch)

    // generate taglist



// Stage 1 - Step 1: BCL to CRAM

    bcl_to_cram(set_up_params.out.bcls, set_up_params.out.lane, set_up_params.out.run_id, set_up_params.out.study_name)//, set_up_params.out.barcodes_file)

// Stage 1 - Step 2: CRAM to BAM
    bcl_to_cram.out
        .flatMap()
        .map {
            basename = it.getBaseName().replaceFirst(/\..*$/, '')  // 1234_5#6.cram -> 1234_5#6
            tag = basename.replaceFirst(/.*#/, '')  // 1234_5#6 -> 6
            [tag, it]
        }.join(all_manifest_data)
        .multiMap {
            tag: it[0]
            cram: it[1]
        }.set { cram_to_bam_input }

    cram_to_bam(cram_to_input.tag, cram_to_input.cram, set_up_params.reference_files, all_manifest_data)
}

// -------------- Check if everything went okay -------------------------------
workflow.onComplete {
    if (workflow.success) {
        log.info """
            ===========================================
            ${ANSI_GREEN}Finished in ${workflow.duration}
            See the report here ==> ${ANSI_RESET}/SOMEDIR/XXX_report.html
            """
            .stripIndent()
    } else {
        log.info """
            ===========================================
            ${ANSI_RED}Finished with errors!${ANSI_RESET}
            """
            .stripIndent()
    }
}
