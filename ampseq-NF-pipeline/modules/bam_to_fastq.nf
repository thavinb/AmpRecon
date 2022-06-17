
process bam_to_fastq {
    /*
    * convert BAM files to FASTQ.
    */
    container ''

    input:
        tuple val(tag), path(input_file)

    output:
        tuple val("${tag}"), path("${base_name}.fastq"), path(${input_file})

    script:
        base_name=input_file.baseName
        """
        bamtofastq \
            < "${base_name}.clipped.bam" \
            > "${base_name}.fastq"
        """
}
