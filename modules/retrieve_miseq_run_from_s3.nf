process retrieve_miseq_run_from_s3 {
    /**
    * Downloads and unzips .tar.gz Miseq run from an S3 bucket.
    */

    input:
        val(uuid_id)

    output:
        tuple val("${params.run_id}"), val("${output_path}"), val("${params.lane}"), val("${params.study_name}"), val("${params.read_group}"), val("${params.library}"), emit: tuple_ch

    script:
        output_path = "${params.results_dir}/${uuid_id}"
        bucket = params.s3_bucket_input

        """
	mkdir -p "${params.results_dir}${uuid_id}"
        cd "${params.results_dir}/${uuid_id}"

        s3cmd get s3://${bucket}/${uuid_id}
        unzip -d ./ "${uuid_id}"
        """
}



