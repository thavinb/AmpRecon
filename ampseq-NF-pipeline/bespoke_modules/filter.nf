process filter {
    /*
    *                                           
    */

    input:
        tuple path(fastq_file)

    output:   
        tuple path("*.filtered.fastq")

    script:
        filtered_fastq_file = "${fastq_file.simpleName}.filtered.fastq"
        """
        mv ${fastq_file} ${filtered_fastq_file}
        """
}
