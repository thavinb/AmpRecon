// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

params.speciation_output_filename = "speciation_output.tsv"

process grc_speciate {
    label "grc_tools"

    input:
        path(genotype_file)
        path(barcodes_file)

    output:
        path("${speciation_out_filename}")

    script:
        grc_settings = params.grc_settings_file_path
        speciation_out_filename = params.speciation_output_filename

        """
        grc_speciate.py \
            --genotype_files ${genotype_file} \
            --output_file "${speciation_out_filename}" \
            --barcodes_file "${barcodes_file}" \
            --config "${grc_settings}"
        """
}
