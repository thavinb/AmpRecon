process unmap {
    /*
    * 
    */

    input:
        path(bam_file)

    output:
        path("*.fastq")

    script:
        fastq_file = "${bam_file.simpleName}.fastq"
        """
        mv ${bam_file} ${fastq_file}
        """
}

