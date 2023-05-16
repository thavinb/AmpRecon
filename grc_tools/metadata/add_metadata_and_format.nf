params.grc1_metadata_filename = params.run_id + "_GRC.txt"
params.grc2_metadata_filename = params.run_id + "_GRC2.txt"
params.barcode_split_formatted_filename = params.run_id + "_Barcodes.txt"

process add_metadata_and_format {
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"

    label "grc_tools"
    input:
        path(manifest_file)
        path(grc1_file)
        path(grc2_file)
        path(barcodes_file)

    output:
        path("${params.grc1_metadata_filename}"), emit: grc1
        path("${params.grc2_metadata_filename}"), emit: grc2
        path("${params.barcode_split_formatted_filename}"), emit: barcodes

    script:
        output_file_name = params.barcode_output_filename

        """
        grc_metadata_and_formatting.py \
            --manifest_file ${manifest_file} \
            --grc1 ${grc1_file} \
            --grc2 ${grc2_file} \
            --barcodes_file ${barcodes_file} \
            --output_file_grc1 "${params.grc1_metadata_filename}" \
            --output_file_grc2 "${params.grc2_metadata_filename}" \
            --output_file_barcodes "${params.barcode_split_formatted_filename}"
        """
}
