process remap {
    /*
    *                                           
    */

    input:
        path(filtered_fastq_file)

    output:   
        path("*.remapped.bam") 

    script:
        remapped_bam_file = "${filtered_fastq_file.simpleName}.remapped.bam"
        """
        mv ${filtered_fastq_file} ${remapped_bam_file}
        """
}
