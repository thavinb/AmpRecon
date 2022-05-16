process cram_to_bam {
    /*
    *                                           
    */

    input:
        path(cram_file)

    output:   
        path("*.bam")

    script:
        bam_file = "${cram_file.simpleName}.bam"
        """
        mv ${cram_file} ${bam_file}
        """
}
