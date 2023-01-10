process bam_to_cram {
    /*
    * split BAM by read group into CRAM.
    */
    errorStrategy 'finish'
    //publishDir "${params.results_dir}/", overwrite: true
    label 'samtools'
    input:
        val(run_id)
        path(adapters_bam_file)
        path(metrics_bam_file)

    output:
        //tuple val(run_id), path("*.cram")
        val(run_id), emit:run_id
        path("*.cram"), emit:cram_fls
        // GAMBIARRA ALERT ----------------------------------------------------
        // this is produced on a step bam_find_adapter step, is here just to
        // simplefy the input for rename_cram_fls =/
        path(metrics_bam_file), emit:metrics_bam_file
        // --------------------------------------------------------------------
    script:
        """
        samtools split \
            --threads 4 \
            -v \
            --output-fmt "cram",no_ref=1 \
            -f "%!.cram" \
            ${adapters_bam_file}
        """
}

