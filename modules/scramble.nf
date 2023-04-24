process scramble_cram_to_bam {
    /**
    * Converts a cram to bam.
    */

    label 'staden'

    input:
        tuple val(file_id), path(cram_file)
    output:
        tuple val(file_id), path("${bam_file}")

    script:
        base_name=cram_file.getSimpleName()
        bam_file="${base_name}.bam"

        """
        scramble \
            -t 1 \
            -7 \
            -I cram \
            -O bam \
            < "${cram_file}" \
            > "${bam_file}"
        """
}

