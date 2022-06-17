params.bam12split_verbose = 0
params.bam12split_level = 0

process bam_split {
    /*
    * Splits BAM rank pairs to single ranks per read
    */
    container ''

    input:
        path(input_bam)

    output:
        tuple val(tag), path("*.split.bam")

    script:
        """
        bam12split verbose=${params.bam12split_verbose} level=${params.bam12split_level} \
            < "${input_bam}" \
            > "${base_name}.split.bam"

        """
}
