#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process make_samplesheet_manifest {

    input:
        val(run_id)
        path(bcl_dir)

    output:
        path("${manifest}")

    script:
        samplesheet = "${bcl_dir}/SampleSheet.csv"
        manifest = "${run_id}_manifest.csv"
        """
        samplesheet_parser.py ${samplesheet} -r ${run_id} -o ${manifest}
        """
}


