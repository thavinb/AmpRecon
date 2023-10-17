// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

params.bambi_compression_level=0

process decode_multiplexed_bam {
    /*
    * Decodes a multiplexed BAM file.
    */

    input:
        path(i2b_output)
        path(tag_file)

    output:
        path("${decoded_bam_file}"), emit: decoded_bam_file
        path("${metrics_file}"), emit: bam_metrics_file

    script:
        compression_level=params.bambi_compression_level
        base_name=i2b_output.baseName
        decoded_bam_file="${base_name}_bambi_decode.bam"
        metrics_file="${decoded_bam_file}.metrics"

        """
        bambi decode \
            --input-fmt="bam" \
            --metrics-file=${metrics_file} \
            --barcode-file=${tag_file} \
            --output=${decoded_bam_file} \
            --output-fmt="bam" \
            --compression-level=${compression_level} \
            ${i2b_output}
        """
}
