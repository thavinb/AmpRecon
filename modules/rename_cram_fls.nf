process rename_cram_fls {
    /*
    * resets a BAM file to a pre-aligned state
    */
    //publishDir "${params.results_dir}", overwrite: true
    label 'pythonBox'
    
    input:
        val(run_id)
        path(bam_metrics)
        path(cram_fls)
        val(lane)
        path(manifest_file)
    output:
        path("*-.cram")

    script:
        // GAMBIARRA ALIERT ---------------------------------------------------
        // this expects the manifest to be at results dir before running this
        // proc
        // ---------------------------------------------------------------------
        """
        renameSamplesCram.py \
                --manifest ${manifest_file} \
                --bam_metrics ${bam_metrics} \
                --cram_file ${cram_fls} \
                --run_id ${run_id} \
                --lane ${lane}
        """
}

