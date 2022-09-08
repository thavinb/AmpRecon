#!/usr/bin/env nextflow

include { download_test_reset_bam_from_s3 } from '../modules/download_test_data.nf'
include { clip_adapters } from '../modules/clip_adapters.nf'
include { download_test_clipped_bam_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	test_bam = download_test_reset_bam_from_s3("21045_1_ref")
	test_sample_id = "123_1"

	under_test = clip_adapters(test_sample_id, test_bam)
	test_bam = under_test.clipped_bam	

	check_md5sum(test_bam, "3f60ca10b6b57d7c4b722b8d9f9b0bbf")
	 	



}

