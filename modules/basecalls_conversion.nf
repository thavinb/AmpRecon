// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

params.bambi_threads = 8
params.bambi_compression_level = 9

process basecalls_conversion {
    /*
    Converts Illumina BCL sequencing run data into a BAM file.
    */

    input:
        val(batch_id)
        path(base_dir)
        val(study_name)

    output:
        path("${output_file}")

    script:
        intensity_dir="${base_dir}/Data/Intensities"
        basecalls_dir="${base_dir}/Data/Intensities/BaseCalls"
        threads="${params.bambi_threads}"
        compression_level="${params.bambi_compression_level}"
        output_file="${batch_id}_bambi_i2b.bam"
        if (!(params.DEBUG_tile_limit == null)){
            tile_limit_arg="--tile-limit=${params.DEBUG_tile_limit}"
        } else {
            tile_limit_arg=""
        }

        """
        bambi i2b \
            --intensity-dir=${intensity_dir} \
            --basecalls-dir=${basecalls_dir} \
            --lane=1 \
            --read-group-id=${batch_id}_1 \
            --platform-unit=${batch_id}_1 \
            --study-name=${study_name} \
            --threads=${threads} \
            --output-file=${output_file} \
            --output-fmt="bam" \
            --compression-level=${compression_level} ${tile_limit_arg}
        """
}
