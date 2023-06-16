process validate_manifest {

    input:
        path(manifest)
        val(panel_names_list)

    script:
        formatted_panel_names_list = panel_names_list.join(" ")
        """
        validate_manifest.py -m ${manifest} -p ${formatted_panel_names_list}
        """
}
