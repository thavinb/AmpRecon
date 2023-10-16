process upload_pipeline_output_to_s3 {
    /*
    Uploads given files to an S3 bucket.
    */
    input:
        path(files)
        val(file_prefix)

    script:
        bucket = params.s3_bucket_output
        """
        s3cmd put "${files}" s3://"${bucket}"/"${params.s3_uuid}"/"${file_prefix}/" --follow-symlinks
        """
}
