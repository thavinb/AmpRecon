process fetch_from_irods {
    /*
    *                                           
    */

    input:
        path(irods_path)

    output:   
        path("*.bam")

    script:
        bam_file = "${irods_path.simpleName}.bam"
        """
        mv ${irods_path} ${bam_file}
        """
}
