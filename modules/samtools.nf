#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

process samtools_sort {
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_tag), path(input_bam)

    output:
        val("${output_directory}"), emit: bam_dir
        tuple val(sample_tag), path("${bam_name}"), emit: bam

    script:
        output_directory = "${params.bam_dir}"
        base_name = input_bam.simpleName
        bam_name="${base_name}.bam"

        """
        samtools sort --threads 2 -o "${bam_name}" "${input_bam}"
        """
}

process samtools_index {
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_tag), path(input_bam)

    output:
        val("${output_directory}"), emit: bam_dir
        tuple val(sample_tag), path(input_bam), path("${bam_name}.bai"), emit: files

    script:
        output_directory = "${params.bam_dir}"
        base_name = input_bam.simpleName
        bam_name="${base_name}.bam"

        """
        samtools index -b "${bam_name}"
        """
}

