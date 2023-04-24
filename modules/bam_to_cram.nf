process bam_to_cram {
    /*
    * split BAM by read group into CRAM.
    */

    label 'samtools'
    input:
        path(adapters_bam_file)

    output:
        path("*.cram")

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

