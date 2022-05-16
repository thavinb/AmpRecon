process fetch_from_irods {
    /*
    *                                           
    */

    input:
        tuple path(irods_path)

    output:   
        tuple path("*.bam")

    script:
        bam_file = "${irods_path.simpleName}.bam"
        """
        mv ${irods_path} ${bam_file}
        """
}
