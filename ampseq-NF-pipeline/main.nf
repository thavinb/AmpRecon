#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
// - workflows
include { prepare_reference } from './pipeline_workflows/step0.1b-prepare-reference.nf'
include { bcl_to_cram } from './pipeline_workflows/step1.1-bcl-to-cram.nf'
include { cram_to_bam } from './pipeline_workflows/step1.2-cram-to-bam.nf'
// - process to extract and validade information expected based on input params
include { get_taglist_file } from './modules/manifest2tag.nf'
include { validate_parameters; load_manifest_ch } from './modules/inputHandling.nf'
include { make_samplesheet_manifest } from './modules/samplesheet_parser.nf'
include { validate_samplesheet_manifest } from './modules/samplesheet_manifest_validation.nf'
include { irods_retrieve } from './modules/irods_retrieve.nf'

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
         --manifest        (required) : ${params.manifest}
         --reference_fasta (required) : ${params.reference_fasta}
         --results_dir                : ${params.results_dir}
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

    // iRODS manifest check and parsing
    if (params.irods_manifest){
        Channel.fromPath(params.irods_manifest, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.file_path))}
            .set{ irods_ch }
    } else {
        Channel.empty().set{ irods_ch }
    }

    // handle reference fasta
    prepare_reference(params.reference_fasta)
    reference_idx_fls = prepare_reference.out

    // Stage 1 - Step 1: BCL to CRAM
    bcl_to_cram(step1_Input_ch)//, set_up_params.out.barcodes_file)
    manifest_step1_1_Out_ch = bcl_to_cram.out.multiMap { it -> run_id: it[0]
                                                                  mnf: it[1]}

    // Stage 1 - Step 2: CRAM to BAM
    cram_to_bam( manifest_step1_1_Out_ch.mnf,
                 reference_idx_fls.bwa_index_fls,
                 reference_idx_fls.fasta_index_fl,
                 reference_idx_fls.dict_fl)

    // Retrieve BAM files from iRODS
    irods_retrieve(irods_ch)

    // Concatenate in-country channel with iRODS channel
    cram_to_bam.out.concat(irods_retrieve.out).set{ bam_files_ch }

    /*
    validate_samplesheet_manifest.out
        .splitCsv(
            header: ["lims_id", "sims_id", "index", "ref", "barcode_sequence", "well", "plate"],
            skip: 18
        ).map { row ->
            ref_map = get_reference_files(row.ref)
            manifest_map = ["lims_id": row.lims_id, "sims_id": row.sims_id, "index": row.index, "ref": row.ref,
                            "barcode_sequence": row.barcode_sequence, "well": row.well, "plate": row.plate]
            all_data_map = manifest_map + ref_map
            indexed_all_data_map = [row.index, all_data_map]
            indexed_all_data_map
        }.set { all_manifest_data }
    */

    /*
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
    */
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
