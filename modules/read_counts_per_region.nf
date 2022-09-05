#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

process read_count_per_region {
    stageInMode 'copy'
    publishDir "${params.results_dir}/", overwrite: true

    input:
        val(run_id)
        path(manifest_file)
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

        grep ${qc_run_id} "${manifest_file}" | awk 'BEGIN {FS=","; OFS=","} {print \$1}' > "${plex_file}"
        python3 ${projectDir}/bin/count_reads_per_region.py \
            --design_file "${qc_cnf_file}" \
            --plex_file "${plex_file}" \
            --input_dir "${bam_directory}" \
            --output "${output_file}"
        """
}

process bam_ref_ch_to_csv {
  input:
    tuple val(bam_name), path(reference_files)
$/
#!/usr/bin/python3
from pathlib import Path

# setup inputs
bam_name = "${bam_name}"
reference_files = "${reference_files}"
publishDir = f"${launchDir}/"

# if manifest already exists, just append new lines
path_to_mnf = f"{publishDir}/bam_ref_ch.csv"
if Path(path_to_mnf).is_file():
    out_mnf = open(f"{path_to_mnf}", "a")

# if manifest does not exist, create file and write header
else:
    out_mnf = open(f"{path_to_mnf}", "w")
    out_mnf.write("sample_tag,reference_fasta\n")

# write manifest line for the bam file
out_mnf.write(f"{bam_name},{reference_files}\n")
out_mnf.close()
/$
}

