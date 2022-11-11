#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

process samtools_sort {
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    input:
        tuple val(sample_tag), path(input_bam)

    output:
        tuple val(sample_tag), path("${bam_name}")

    script:
        base_name = input_bam.baseName
        bam_name="${base_name}.bam"

        """
        echo ${input_bam}
        samtools sort --threads 2 -o "${bam_name}" "${input_bam}"
        """
}

process samtools_index {
    // saveAs argument to ensure that BQSR saved output has same name as when BQSR not run
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    input:
        tuple val(sample_tag), path(input_bam)

    output:
        tuple val(sample_tag), path(input_bam), path("${bam_name}.bai")

    script:
        base_name = input_bam.baseName
        bam_name="${base_name}.bam"

        """
        samtools index -b "${bam_name}"
        """
}
