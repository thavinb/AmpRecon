process unmap {
    /*
    * 
    */

    input:
        tuple path(bam_file)

    output:
        tuple path("*.fastq")

    script:
        fastq_file = "${bam_file.simpleName}.fastq"
        """
        mv ${bam_file} ${fastq_file}
        """
}

