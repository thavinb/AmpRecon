process validate_samplesheet_manifest {
    
    input:
        tuple val(run_id), path(manifest)
        val(panel_names_list)
    output:
        tuple val(run_id), path("${manifest}")

    script:
        formatted_panel_names_list = "${panel_names_list}".replaceAll("[\\,\\[\\]]", "")
        """
        manifest_validation.py -m ${manifest} -p ${formatted_panel_names_list}
        """
    }
