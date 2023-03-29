// enable dsl2
nextflow.enable.dsl = 2

process read_count_per_region {
    stageInMode 'copy'
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    label 'python_plus_samtools'

    input:
        path(bam_file_list)
        path(bam_files_and_indices)
        tuple val(panel_name), file(annotation_file)

    output:
        path("${output_file}"), emit: qc_csv_file

    script:
        output_file = "${panel_name}_reads_per_region.csv"
        plex_file = "${panel_name}.plex"

        """
        grep ${panel_name} "${bam_file_list}" | awk 'BEGIN {FS=","; OFS=","} {print \$1}' > "${plex_file}"
        count_reads_per_region.py  \
            --design_file "${annotation_file}" \
            --plex_file "${plex_file}" \
            --input_dir "." \
            --output "${output_file}"
        """
}

process files_and_panels_to_csv {
  label "pythonBox"
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

