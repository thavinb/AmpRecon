// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

params.barcode_output_name = "barcode.intermediate.tsv"
params.barcode_intermediate_splitout_name = "barcoding_output.split_out.intermediate.tsv"

process grc_barcoding {
    label "py_pandas"
    label "grc_tools"
    input:
        path(genotype_file)

    output:
        path("${output_file_name}"), emit: barcoding_file
        path("${split_output_file_name}"), emit: barcoding_split_out_file
    script:
        grc_settings = params.grc_settings_file_path
        output_file_name = params.barcode_output_name
        split_output_file_name = params.barcode_intermediate_splitout_name

        """
        grc_barcoding.py \
            --genotype_files ${genotype_file} \
            --output_file "${output_file_name}" \
            --config "${grc_settings}" \
            --output_file_split_out ${split_output_file_name}
        """
}
