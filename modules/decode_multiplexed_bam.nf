params.bambi_compression_level=0

process decode_multiplexed_bam {
    /*
    * Decodes a multiplexed BAM file.
    */
    //publishDir "${params.results_dir}/", overwrite: true
    
    input:
        tuple val(run_id), path(base_dir), val(lane), val(study_name), val(read_group), val(library), path(barcode_file), path(i2b_output)
        //path(input_file) // ex. bambi_i2b.bam
        //path(barcode_file) // ex.

    output:
        tuple val(run_id), path("${decoded_bam_file}"), path("${metrics_file}")
        //path("${output_file}"), emit: output_file
        //path("${metrics_file}"), emit: metrics_file

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
