params.scramble = 'scramble'

process scramble_sam_to_bam {
    /**
    * Converts a sam to bam.
    */

    

    input:
        val(tag)
        path(sam_file)

    output:
        tuple val(tag), path("${bam_file}")

    script:
        base_name=sam_file.getBaseName()
        bam_file="${base_name}.bam"

        scramble=params.scramble

        """
        ${scramble} -0 -I sam -O bam "${sam_file}" "${bam_file}"
        """
}


params.scramble_cram_to_bam_threads = 11
params.scramble_cram_to_bam_compression_level = -7

process scramble_cram_to_bam {
    /**
    * Converts a cram to bam.
    */

    

    input:
        tuple val(tag), path(cram_file)
    output:
        tuple val(tag), path("${bam_file}")

    script:
        base_name=cram_file.getSimpleName()
        bam_file="${base_name}.bam"

        scramble=params.scramble

        """
        ${scramble} \
            -t ${params.scramble_cram_to_bam_threads} \
            ${params.scramble_cram_to_bam_compression_level} \
            -I cram \
            -O bam \
            < "${cram_file}" \
            > "${bam_file}"
        """
}

