#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_samplesheet_manifest {

    input:
        path(samplesheet_csv) // usually "${bcl_dir}/SampleSheet.csv" 

    output:
        path("${manifest}")

    script:
        manifest = "${params.run_id}_manifest.csv"
        """
        samplesheet_parser.py ${samplesheet_csv} -o ${manifest}
        """
}


