process bcl_to_bam {
    /*
    *                                           
    */

    input:
        tuple path(bcl_file)

    output:   
        tuple path("*.bam")

    script:
        bam_file = "${bcl_file.simpleName}.bam"
        """
        mv ${bcl_file} ${bam_file}
        """
}
