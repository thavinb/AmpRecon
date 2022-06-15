process decode_multiplexed_bam {
    /*
    * Decodes a multiplexed BAM file.
    */
    container ''

    input:
        path(input_file) // ex. bambi_i2b.bam
        path(barcode_file) // ex.

    output:
        path("${output_file}"), emit: output_file
        path("${metrics_file}"), emit: metrics_file

    script:
        compression_level=params.bambi_compression_level
        base_name=${input_file}.baseName
        output_file="${base_name}.bambi_decode.bam"
        metrics_file="${output_file}.metrics"

        """
        bambi decode \
            --input-fmt="bam" \
            --metrics-file=${metrics_file} \
            --barcode-file=${barcode_file} \
            --output=${output_file} \
            --output-fmt="bam" \
            --compression-level=${compression_level} \
            ${input_file}
        """
}
