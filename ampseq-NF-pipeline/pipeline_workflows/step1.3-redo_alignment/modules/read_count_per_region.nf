#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

process make_file_list {

    input:
        path(bam_directory)

    output:
        path("file_list.csv")

    script:

        """
        set -eo pipefail

        for file in "${bam_directory}"/*.bam;
        do
            base_name=\$(basename "\$file" .bam)

            echo "\$base_name,PFA_Spec" >> "file_list.csv"
            echo "\$base_name,PFA_GRC1_v1.0" >> "file_list.csv"
            echo "\$base_name,PFA_GRC2_v1.0" >> "file_list.csv"
        done
        """
}

process read_count_per_region {
    stageInMode 'copy'
    publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        val(run_id)
        path(file_list)
        path(bam_directory)
        val(qc_run_id)
        path(qc_cnf_file)

    output:
        path("${output_file}"), emit: qc_csv_file
        path("${plex_file}"), emit: qc_plex_file

    script:
        output_file = "${run_id}_${qc_run_id}_reads_per_region.csv"
        plex_file = "${run_id}_${qc_run_id}.plex"

        """
        set -eo pipefail

        grep ${qc_run_id} "${file_list}" | awk 'BEGIN {FS=","; OFS=","} {print \$1}' > "${plex_file}"
        python3 ${projectDir}/pipeline_workflows/step1.3-redo_alignment/modules/count_reads_per_region.py \
            --design_file "${qc_cnf_file}" \
            --plex_file "${plex_file}" \
            --input_dir "${bam_directory}" \
            --output "${output_file}"
        """
}
