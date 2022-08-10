#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

process sort_and_index {
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: true
    input:
        tuple val(sample_tag), path(input_bam)

    output:
        val("${output_directory}"), emit: bam_dir
        path("${bam_name}"), emit: bam

    script:
        output_directory = "${params.bam_dir}"
        base_name = aligned_sam.simpleName
        bam_name="${base_name}.sorted.bam"

        """
        samtools sort --threads 2 -o "${bam_name}" "${input_bam}"
        samtools index -b "${bam_name}"
        """
}

process read_count_per_region_qc {
    stageInMode 'copy'

    input:
        val(key)
        path(bam_directory)
        val(qc_run_id)
        path(qc_cnf_file)

    output:
        path("${output_file}"), emit: qc_csv_file
        path("${plex_file}"), emit: qc_plex_file

    script:
        output_file = "${key}_${qc_run_id}_reads_per_region.csv"
        plex_file = "${key}_${qc_run_id}.plex"

        """
        set -eo pipefail

        for file in "${bam_directory}"/*.bam;
        do
            base_name=\$(basename "\$file" .bam)

            echo "\$base_name,PFA_Spec" >> "mock_manifest.csv"
            echo "\$base_name,PFA_GRC1_v1.0" >> "mock_manifest.csv"
            echo "\$base_name,PFA_GRC2_v1.0" >> "mock_manifest.csv"
        done

        grep ${qc_run_id} "mock_manifest.csv" | awk 'BEGIN {FS=","; OFS=","} {print \$1}' > "${plex_file}"
        python3 ${projectDir}/pipeline_workflows/step1.3-redo_alignment/modules/count_reads_per_region.py \
            --design_file "${qc_cnf_file}" \
            --plex_file "${plex_file}" \
            --input_dir "${bam_directory}" \
            --output "${output_file}"
        """
}
