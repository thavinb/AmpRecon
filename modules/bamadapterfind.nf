process BAMADAPTERFIND {
    /*
    * Searches for sequencing adapter contamination with a BAM  file.
    */
    publishDir "${params.results_dir}/${run_id}/", overwrite: true

    input:
        tuple val(run_id), path(bam_file)

    output:
        //tuple val(run_id), path("${adapters_bam_file}")
        tuple val(run_id), path("${adapters_bam_file}"), path("${metrics_bam_file}")

    script:
        adapters_bam_file = "${bam_file}.adapters"
        metrics_bam_file = "${bam_file}.metrics"
        // TODO add threads option ?
        """
        bamadapterfind level=9 < ${bam_file} > ${adapters_bam_file}
        """

}
