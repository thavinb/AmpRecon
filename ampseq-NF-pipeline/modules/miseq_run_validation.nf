#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process miseq_run_validation {
    /*
    * Checks the following exist for a given MiSeq run: SampleSheet.csv, `Basecalls/L001/ and .bcl files.
    */

    publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

    input:
        tuple val(run_id), val(bcl_dir), val(not_used1), val(not_used2), val(not_used3), val(not_used4), val(not_used5)

    script:
        samplesheet = "${bcl_dir}/SampleSheet.csv"
        basecalls = "${bcl_dir}/Data/Intensities/BaseCalls/L001/"
        BCLs = file("${bcl_dir}/Data/Intensities/BaseCalls/L001/*/*.bcl")
        BCL_PATH=BCLs.first()

        if (!file("${samplesheet}").exists()) {
            log.error("The samplesheet for run:${run_id} is missing from ${bcl_dir}/.")
        }

        if (!file("${basecalls}").exists()) {
            log.error("Run:${run_id} is missing a BaseCalls/L001/ directory.")
        }

        if (!file(BCL_PATH).exists()) {
            log.error("Run:${run_id} is missing BCL data.")
        }
        """
        """
}
