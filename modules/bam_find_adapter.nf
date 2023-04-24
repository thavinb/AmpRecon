process bam_find_adapter {
    /*
    * Searches for sequencing adapter contamination with a BAM  file.
    */
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    
    input:
        tuple val(run_id), path(bam_file), path(not_used1)

    output:
        //tuple val(run_id), path("${adapters_bam_file}")
        val(run_id), emit: run_id
        path("${adapters_bam_file}"), emit: bam_adapter_file
        path("${metrics_bam_file}"), emit: bam_metrics_file

    script:
        adapters_bam_file = "${bam_file}.adapters"
        metrics_bam_file = "${bam_file}.metrics"
        // TODO add threads option ?
        """
        bamadapterfind level=9 < ${bam_file} > ${adapters_bam_file}
        """

}

