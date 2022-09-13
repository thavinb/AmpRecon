// enable dsl2
nextflow.enable.dsl = 2

process read_count_per_region {
    stageInMode 'copy'
    publishDir "${params.results_dir}/", overwrite: true

    input:
        val(run_id)
        path(bam_file_list)
        path(bam_files_and_indices)
        tuple val(pannel_name), file(annotation_file)
        //val(qc_run_id)
        //path(qc_cnf_file)

    output:
        path("${output_file}"), emit: qc_csv_file
        path("${plex_file}"), emit: qc_plex_file

    script:
        output_file = "${run_id}_${pannel_name}_reads_per_region.csv"
        plex_file = "${run_id}_${pannel_name}.plex"

        """
        grep ${pannel_name} "${bam_file_list}" | awk 'BEGIN {FS=","; OFS=","} {print \$1}' > "${plex_file}"
<<<<<<< HEAD
        python3 ${projectDir}/pipeline_workflows/step1.3-redo_alignment/modules/count_reads_per_region.py \
=======
        python3 ${projectDir}/modules/count_reads_per_region.py \
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
            --design_file "${annotation_file}" \
            --plex_file "${plex_file}" \
            --input_dir "." \
            --output "${output_file}"
        """
}

process files_and_panels_to_csv {
<<<<<<< HEAD
  input:
    val(file_names_list)
  output:
    path("file_names_panel_list.csv")
=======
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


process bam_ref_ch_to_csv {
  input:
    tuple val(sample_tag), path(reference_files)
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
$/
#!/usr/bin/python3
from pathlib import Path

<<<<<<< HEAD
path_to_mnf = "file_names_panel_list.csv"
names_list = list("${file_names_list}".strip("[]").replace(" ", "").split(","))
=======
# setup inputs
sample_tag = "${sample_tag}"
reference_files = "${reference_files}"
publishDir = f"${launchDir}/"
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e

out_mnf = open(f"{path_to_mnf}", "w")
out_mnf.write("file_name\n")

for file_name in names_list:

<<<<<<< HEAD
    out_mnf.write(f"{file_name}\n")
=======
# write manifest line for the bam file
out_mnf.write(f"{sample_tag},{reference_files}\n")
>>>>>>> ecf44664ae9377f7ff4a34dcb4206db2d38aaf9e
out_mnf.close()
/$
}

