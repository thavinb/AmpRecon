#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include { reads_retrieval_workflow  } from './pipeline_workflows/reads_retrieval_workflow.nf'
include { filtering_workflow } from './pipeline_workflows/filtering_workflow.nf'
include { variant_calling_workflow } from './pipeline_workflows/variant_calling_workflow.nf'
include { genotyping_workflow } from './pipeline_workflows/genotyping_workflow.nf'

workflow set_up_params {
    main:
        Channel.fromPath(params.bcl_path).set {bcls}
        Channel.fromPath(params.irods_paths).set {irods_paths}
    emit:
        bcls
        irods_paths
}

// Main entry-point workflow

workflow {

    set_up_params()

// Step 1 - Reads retrieval

    reads_retrieval_workflow(set_up_params.out)

// Step 2 - Added filtering

    filtering_workflow(reads_retrieval_workflow.out)

// Step 3 - Variant calling

    variant_calling_workflow(filtering_workflow.out)

// Step 4 - Drug resistance genotyping

    genotyping_workflow(variant_calling_workflow.out)
}
