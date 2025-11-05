// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

params.grc_intermediate_name = "GRC.intermediate.tsv"

process grc_assemble {
    label "py_pandas"
    label "grc_tools"
    input:
        path(grc_components_list)

    output:
        path("${output_grc}")

    script:
        output_grc = params.grc_intermediate_name
        grc_component_files = grc_components_list.join(" ")
        """
        grc_assemble.py -grcs_in ${grc_component_files} -grc_out_name ${output_grc}
        """
}
