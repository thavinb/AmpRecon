params.bambi_select_compression_level = 0
params.bambi = "bambi"
process alignment_filter {
    label 'bambi'
    input:
        val(sample_tag)
        path(merged_bam) // e.g. bam

    output:
        val(sample_tag), emit: sample_tag
        path("${selected_bam}"), emit: selected_bam
        path("${output_metrics_file}"), emit: selected_bam_metrics

    script:
        bambi=params.bambi
        base_name=merged_bam.getSimpleName()
        selected_bam="${base_name}.selected.bam"
        output_metrics_file="${base_name}.metrics"

        """
        ${bambi} select \
            --compression-level=${params.bambi_select_compression_level} \
            --input "${merged_bam}" \
            --output "${selected_bam}" \
            -m "${output_metrics_file}"
        """
}
