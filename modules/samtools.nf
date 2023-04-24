process samtools_index {

    input:
        tuple val(sample_tag), path(input_bam)

    output:
        tuple val(sample_tag), path(input_bam), path("${bam_name}.bai")

    script:
        base_name = input_bam.baseName
        bam_name="${base_name}.bam"

        """
        samtools index -b "${bam_name}"
        """
}
