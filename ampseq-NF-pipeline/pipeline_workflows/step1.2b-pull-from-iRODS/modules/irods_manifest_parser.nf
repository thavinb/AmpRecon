process irods_manifest_parser {
    /*
    Creates iRODS file path from amplicon lanes file and retrieves sample ID from iRODS.
    */
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true

    input:
         tuple val(id_run), val(WG_lane)

    output:
        tuple env(sample_id), val("${iRODS_file_path}"), val(id_run)

    script:
        iRODS_file_path = "/seq/${id_run}/${WG_lane}.cram"
        """
        imeta ls -d "${iRODS_file_path}" > imeta_data
        sample_id=\$(cat imeta_data | grep -A1 "attribute: sample_id" | tail -n 1 | cut -d" " -f2)
        """

}

