#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_samplesheet_manifest {
    label 'pythonBox'
    input:
        tuple val(run_id), path(bcl_dir)

    output:
        tuple val(run_id), path("${manifest}"), emit: tuple
        path("${manifest}"), emit: manifest_file

    script:
        samplesheet = "${bcl_dir}/SampleSheet.csv"
        manifest = "${run_id}_manifest.csv"
        """
        samplesheet_parser.py ${samplesheet} -r ${run_id} -o ${manifest}
        """
}


