process bcl_to_bam {
    /*
    *                                           
    */

    input:
        path(bcl_file)

    output:   
        path("*.bam")

    script:
        bam_file = "${bcl_file.simpleName}.bam"
        """
        mv ${bcl_file} ${bam_file}
        """
}
