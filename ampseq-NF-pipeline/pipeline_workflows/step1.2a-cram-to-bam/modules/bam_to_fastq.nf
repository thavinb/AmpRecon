
process bam_to_fastq {
    /*
    * convert BAM files to FASTQ.
    */
    //publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        tuple val(sample_tag), path(clipped_bam)

    output:
        tuple val(sample_tag), path("${base_name}.fastq")

    script:
        base_name=clipped_bam.baseName
        """
        bamtofastq \
            < "${clipped_bam}" \
            > "${base_name}.fastq"
        """
}
