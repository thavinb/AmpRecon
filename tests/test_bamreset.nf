#!/usr/bin/env nextflow


include { download_test_collated_bam_from_s3 } from '../modules/download_test_data.nf'
include { bam_reset } from '../modules/bam_reset.nf'
include { download_test_reset_bam_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	test_ch = download_test_collated_bam_from_s3("21045_1_ref")

	test_input_ch = generate_test_channel(test_ch)

	under_test = bam_reset(test_input_ch)
	test_bam = under_test.reset_bam	

	check_md5sum(test_bam, "4ebf9a0f0330a59c826ccc9b6e1206f4")

	 	



}



process generate_test_channel {

	input:
	path(bam_file)

	output:
	val("123_1")
        path(bam_file)

	script:
	"""
	echo
	"""

}
