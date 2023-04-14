params.barcode_output_filename = "barcode_output.tsv"

process grc_barcoding {
    label "grc_tools"
    input:
        path(genotype_file)

    output:
        path("${output_file_name}")

    script:
        grc_settings = params.grc_settings_file_path
        output_file_name = params.barcode_output_filename

        """
        grc_barcoding.py \
            --genotype_files ${genotype_file} \
            --output_file "${output_file_name}" \
            --config "${grc_settings}"
        """
}
