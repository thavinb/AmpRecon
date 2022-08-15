process rename_cram_fls {
    /*
    * resets a BAM file to a pre-aligned state
    */
    publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        //path(manifest_csv)
        val(run_id)
        path(bam_metrics)
        path(cram_fls)

    output:
        tuple val(run_id), path("*_.cram")

    script:
        // GAMBIARRA ALIERT ---------------------------------------------------
        // this expects the manifest to be at results dir before running this
        // proc
        manifest_csv = "${params.results_dir}/${run_id}/${run_id}_manifest.csv"
        // ---------------------------------------------------------------------
        """
        python3 ${projectDir}/pipeline_workflows/step1.1-bcl-to-cram/modules/scripts/renameSamplesCram.py \
                --manifest ${manifest_csv} \
                --bam_metrics ${bam_metrics} \
                --cram_file ${cram_fls}
        """
}
