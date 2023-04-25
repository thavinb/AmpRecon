process retrieve_miseq_run_from_s3 {
    /**
    * Downloads and unzips .tar.gz Miseq run from an S3 bucket.
    */

    input:
        val(uuid_id)

    output:
        val("${output_path}")

    when:
        params.s3_uuid != "-1"

    script:
        output_path = "${params.results_dir}/${uuid_id}"
        bucket = params.s3_bucket_input

        """
	    mkdir -p "${params.results_dir}${uuid_id}"
        cd "${params.results_dir}/${uuid_id}"

        s3cmd get --force s3://${bucket}/${uuid_id}
        unzip -o -q -d ./ "${uuid_id}"
        """
}



