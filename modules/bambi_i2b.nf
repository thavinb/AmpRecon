params.bambi_threads = 8
params.bambi_compression_level = 9

process BAMBI_I2B {
    /*
    * Converts Illumina BCL sequencing run data into a BAM file.
    */
    container ''

    input:
        tuple val(run_id), path(base_dir), val(lane), val(study_name), val(read_group)

    output:
        tuple val(run_id), path("${output_file}")

    script:
        intensity_dir="${base_dir}/Data/Intensities"
        basecalls_dir="${base_dir}/Data/Intensities/BaseCalls"
        read_group_id="${read_group}_${lane}"
        threads="${params.bambi_threads}"
        compression_level="${params.bambi_compression_level}"
        output_file="${run_id}_bambi_i2b.bam"

        """
        bambi i2b \
            --intensity-dir=${intensity_dir} \
            --basecalls-dir=${basecalls_dir} \
            --lane=${lane} \
            --read-group-id=${read_group_id} \
            --study-name=${study_name} \
            --threads=${threads} \
            --output-file=${output_file} \
            --output-fmt="bam" \
            --compression-level=${compression_level}
        """
}
