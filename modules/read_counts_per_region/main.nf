// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process READ_COUNT_PER_REGION {
    tag "${panel_name}"
    label 'pysam'
    publishDir "${params.results_dir}/read_counts/", overwrite: true, mode: "copy"

    input:
        tuple val(panel_name), path(dsgn), path(fai), path(bams), path(bais)

    output:
        path("*_reads_per_region.csv"), emit: qc_csv_file

    script:
        def output_file = "${panel_name}_reads_per_region.csv"

        """
        set -e
        set -o pipefail

        find -name '*.bam' | sed 's/.bam//g' > ${panel_name}.plex

        count_reads_per_region.py  \
            --design_file "${dsgn}" \
            --plex_file "${panel_name}.plex" \
            --input_dir "." \
            --output "${output_file}"
        """
}

