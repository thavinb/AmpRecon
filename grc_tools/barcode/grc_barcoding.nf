params.barcode_output_filename = "barcode_output.tsv"

process grc_barcoding {
    label "grc_tools"
    input:
        val(genotype_file_list)

    output:
        path("${output_file_name}")

    script:
        grc_settings = params.grc_settings_file_path
        output_file_name = params.barcode_output_filename
        genotype_file_list = genotype_file_list.join(" ")

        """
        grc_barcoding.py \
            --genotype_files ${genotype_file_list} \
            --output_file "${output_file_name}" \
            --config "${grc_settings}"
        """
}
