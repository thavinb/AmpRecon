process irods_manifest_parser {
    /*
    Creates iRODS file path from amplicon lanes file and retrieves sample ID from iRODS.
    */
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true

    input:
         tuple val(id_run), val(WG_lane), val(irods_path)

    output:
        tuple env(sample_id), val(irods_path), val(WG_lane)

    script:
        //iRODS_file_path = "/seq/${id_run}/${WG_lane}.cram"
        iRODS_file_path = "${irods_path}"
        flnm = "${irods_path}".split('/')[-1]
        WG_lane = "${flnm}".split('\\.')[0] // split does not work with dot (it recognize it as a regular expression)
        //"${irods_path}".split('/')[-1].split('.')[0]
        
        """
        imeta ls -d "${iRODS_file_path}" > imeta_data
        sample_id=\$(cat imeta_data | grep -A1 "attribute: sample_id" | tail -n 1 | cut -d" " -f2)
        """

}
