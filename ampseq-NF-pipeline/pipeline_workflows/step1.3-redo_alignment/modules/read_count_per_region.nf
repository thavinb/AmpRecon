#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

process read_count_per_region {
    stageInMode 'copy'
    publishDir "${params.results_dir}/", overwrite: true

    input:
        val(run_id)
        path(bam_file_list)
        path(bam_files_and_indices)
        val(qc_run_id)
        path(qc_cnf_file)

    output:
        path("${output_file}"), emit: qc_csv_file
        path("${plex_file}"), emit: qc_plex_file

    script:
        output_file = "${run_id}_${qc_run_id}_reads_per_region.csv"
        plex_file = "${run_id}_${qc_run_id}.plex"

        """

        grep ${qc_run_id} "${bam_file_list}" | awk 'BEGIN {FS=","; OFS=","} {print \$1}' > "${plex_file}"
        python3 ${projectDir}/pipeline_workflows/step1.3-redo_alignment/modules/count_reads_per_region.py \
            --design_file "${qc_cnf_file}" \
            --plex_file "${plex_file}" \
            --input_dir "." \
            --output "${output_file}"
        """
}

process files_and_panels_to_csv {
  input:
    val(file_names_list)
  output:
    path("file_names_panel_list.csv")
$/
#!/usr/bin/python3
from pathlib import Path

path_to_mnf = "file_names_panel_list.csv"
names_list = list("${file_names_list}".strip("[]").replace(" ", "").split(","))

out_mnf = open(f"{path_to_mnf}", "w")
out_mnf.write("file_name\n")

for file_name in names_list:

    out_mnf.write(f"{file_name}\n")
out_mnf.close()
/$
}
