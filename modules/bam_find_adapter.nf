process bam_find_adapter {
    /*
    * Searches for sequencing adapter contamination with a BAM  file.
    */
    publishDir "${params.results_dir}/${run_id}/", overwrite: true

    input:
        tuple val(run_id), path(bam_file), path(not_used1)

    output:
        tuple val(run_id), path("${adapters_bam_file}")

    script:
        adapters_bam_file = "${bam_file}.adapters"
        // TODO add threads option ?
        """
        bamadapterfind level=9 < ${bam_file} > ${adapters_bam_file}
        """

}
