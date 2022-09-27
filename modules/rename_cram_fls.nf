process rename_cram_fls {
    /*
    * resets a BAM file to a pre-aligned state
    */
    //publishDir "${params.results_dir}", overwrite: true
    label 'pythonBox'
    
    input:
        //path(manifest_csv)
        val(run_id)
        path(bam_metrics)
        path(cram_fls)
        val(lane)

    output:
        path("*-.cram")

    script:
        // GAMBIARRA ALIERT ---------------------------------------------------
        // this expects the manifest to be at results dir before running this
        // proc
        manifest_csv = "${params.results_dir}/${run_id}_manifest.csv"
        // ---------------------------------------------------------------------
        """
        renameSamplesCram.py \
                --manifest ${manifest_csv} \
                --bam_metrics ${bam_metrics} \
                --cram_file ${cram_fls} \
                --run_id ${run_id} \
                --lane ${lane}
        """
}

