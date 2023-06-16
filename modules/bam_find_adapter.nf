process bam_find_adapter {
    /*
    * Searches for sequencing adapter contamination with a BAM  file.
    */

    input:
        path(bam_file)

    output:
        path("${adapters_bam_file}")

    script:
        adapters_bam_file = "${bam_file}.adapters"
        """
        bamadapterfind level=9 < ${bam_file} > ${adapters_bam_file}
        """

}

