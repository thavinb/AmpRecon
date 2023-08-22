params.grc1_name = params.run_id + "_GRC.txt"
params.grc2_name = params.run_id + "_GRC2.txt"
params.barcodes_name = params.run_id + "_Barcodes.txt"

process add_metadata_and_format {
    publishDir "${params.results_dir}/grcs_barcodes/", overwrite: true, mode: "copy"

    label "grc_tools"
    input:
        path(manifest_file)
        path(grc1_file)
        path(grc2_file)
        path(barcodes_file)

    output:
        path("${output_grc1}"), emit: grc1
        path("${output_grc2}"), emit: grc2
        path("${output_barcodes}"), emit: barcodes

    script:
        output_grc1 = params.grc1_name
        output_grc2 = params.grc2_name
        output_barcodes = params.barcodes_name
        grc_settings = params.grc_settings_file_path

        """
        grc_metadata_and_formatting.py \
            --manifest_file ${manifest_file} \
            --grc1 ${grc1_file} \
            --grc2 ${grc2_file} \
            --barcodes_file ${barcodes_file} \
            --output_file_grc1 "${output_grc1}" \
            --output_file_grc2 "${output_grc2}" \
            --output_file_barcodes "${output_barcodes}" \
            --config ${grc_settings}
        """
}
