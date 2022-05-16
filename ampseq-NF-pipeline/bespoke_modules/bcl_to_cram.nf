process bcl_to_cram {
    /*
    *
    */

    input:
        path(bcl_file)

    output:
        path("*.cram")

    script:
        cram_file = "${bcl_file.simpleName}.cram"
        """
        mv ${bcl_file} ${cram_file}
        """
}
