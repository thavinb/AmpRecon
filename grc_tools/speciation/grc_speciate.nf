params.speciation_output_filename = "speciation_output.tsv"

process grc_speciate {
    label "grc_tools"

    input:
        val(genotype_file_list)
        path(barcodes_file)

    output:
        path("${speciation_out_filename}")

    script:
        grc_settings = params.grc_settings_file_path
        speciation_out_filename = params.speciation_output_filename
        genotype_file_list = genotype_file_list.join(" ")

        """
        grc_speciate.py \
            --genotype_files ${genotype_file_list} \
            --output_file "${speciation_out_filename}" \
            --barcodes_file "${barcodes_file}" \
            --config "${grc_settings}"
        """
}
