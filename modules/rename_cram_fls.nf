process rename_cram_fls {
    /*
    * resets a BAM file to a pre-aligned state
    */
    
    input:
        val(run_id)
        path(bam_metrics)
        path(cram_fls)
        val(lane)
        path(manifest_file)

    output:
        path("*_*#*_*.cram")

    script:
        """
        renameSamplesCram.py \
                --manifest ${manifest_file} \
                --bam_metrics ${bam_metrics} \
                --cram_file ${cram_fls} \
                --run_id ${run_id} \
                --lane ${lane}
        """
}

