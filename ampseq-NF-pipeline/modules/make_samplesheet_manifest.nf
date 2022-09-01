#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_samplesheet_manifest {
    publishDir "${params.results_dir}/", mode: 'copy', overwrite: true

    input:
    // TODO: only get values that are used
    //path samplesheet
    tuple val(run_id), path(bcl_dir), val(not_used1), val(not_used2), val(not_used3), val(not_used4), val(not_used5)
    //val(run_id)
    //val(bcl_dir)
    output:
    tuple val(run_id), path("${run_id}_manifest.csv"), emit: tuple
    path("${run_id}_manifest.csv"), emit: manifest_file
    script:
    samplesheet = "${bcl_dir}/SampleSheet.csv"
    manifest = "${bcl_dir}/${run_id}_manifest.csv"
    """
    python3 ${workflow.projectDir}/modules/samplesheet_parser.py ${samplesheet} -r ${run_id}
    ln ${manifest} ./
    """
}
