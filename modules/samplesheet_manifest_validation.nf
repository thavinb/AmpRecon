process validate_samplesheet_manifest {
    label 'pythonBox'
    input:
        path(manifest)
        val(panel_names_list)

    script:
        formatted_panel_names_list = "${panel_names_list}".replaceAll("[\\,\\[\\]]", "")
        """
        manifest_validation.py -m ${manifest} -p ${formatted_panel_names_list}
        """
    }
