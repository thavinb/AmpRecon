process bam_to_cram {
    /*
    * split BAM by read group into CRAM.
    */
    container ''

    input:
        path(input_bam)

    output:
        path("*.cram")

    script:
        """
        samtools split \
            --threads 4 \
            -v \
            --output-fmt "cram",no_ref=1 \
            -f "%!.cram" \
            ${input_bam}
        """
}
