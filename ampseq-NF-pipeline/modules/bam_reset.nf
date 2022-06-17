params.bamreset_resetaux = 0
params.bamreset_level = 9
params.bamreset_verbose = 0


process bam_reset {
    /*
    * resets a BAM file to a pre-aligned state
    */
    container ''

    input:
        tuple val(tag), path(input_file)

    output:
        tuple val("${tag}"), path("${base_name}.prealigned.bam")

    script:
        base_name=input_file.baseName
        """
        bamreset \
            resetaux=${params.bamreset_resetaux} \
            level=${params.bamreset_level} \
            verbose=${params.bamreset_verbose} \
            < "${input_file}" \
            > "${base_name}.prealigned.bam"

        """
}
