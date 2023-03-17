params.grc_1_name = "GRC1.tsv"

process grc_assemble {
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    label "grc_tools"
    input:
        val(grc_components_list)

    output:
        path("${output_grc_1}")

    script:
        output_grc_1 = params.grc_1_name
        grc_components_list = grc_components_list.join(" ")
        """
        grc_assemble.py -grcs_in ${grc_components_list} -grc_out_name ${output_grc_1}
        """
}
