params.bambi_compression_level=0

process BAMBI_DECODE {
    /*
    * Decodes a multiplexed BAM file.
    */
    publishDir "${params.results_dir}/${run_id}/", overwrite: true

    input:
        tuple val(run_id), path(i2b_output)
        path(barcode_file)

    output:
        tuple val(run_id), path("${decoded_bam_file}"), path("${metrics_file}")

    script:
        compression_level=params.bambi_compression_level
        base_name=i2b_output.baseName
        decoded_bam_file="${base_name}_bambi_decode.bam"
        metrics_file="${decoded_bam_file}.metrics"

        """
        bambi decode \
            --input-fmt="bam" \
            --metrics-file=${metrics_file} \
            --barcode-file=${barcode_file} \
            --output=${decoded_bam_file} \
            --output-fmt="bam" \
            --compression-level=${compression_level} \
            ${i2b_output}
        """
}
