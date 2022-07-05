params.bambi_select_compression_level = 0
params.bambi = "bambi"

process bambi_select {
    /**
    *
    */
    input:
        val(sample_tag)
        path(input_file) // e.g. bam

    output:
        tuple val(sample_tag), path("${output_file}"), path("${output_metrics_file}")
        //val(sample_tag), emit: sample_tag
        //path("${output_file}"), emit: selected_bam
        //path("${output_metrics_file}"), emit: selected_bam_metrics

    script:
        bambi=params.bambi
        base_name=input_file.getBaseName()
        output_file="${base_name}.selected.bam"
        output_metrics_file="${base_name}.metrics"

        """
        ${bambi} select \
            --compression-level=${params.bambi_select_compression_level} \
            --input "${input_file}" \
            --output "${output_file}" \
            -m "${output_metrics_file}"
        """
}
