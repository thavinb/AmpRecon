#!/usr/bin/env nextflow

include { download_test_aligned_sam_from_s3 } from '../modules/download_test_data.nf'
include { scramble_sam_to_bam } from '../modules/scramble.nf'
include { download_test_scrambled_bam_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	test_ch = download_test_aligned_sam_from_s3("21045_1_ref")
	test_tag = "123"

	under_test = scramble_sam_to_bam(test_tag, test_ch)
	test_bam = under_test.map { it -> it[1] }
	
	check_md5sum(test_bam, "a35101ba730765b00d6a8bb426dd3576")


}



process generate_test_channel {

	input:
	path(sam_file)

	output:
	val("123_1")
        path(sam_file)

	script:
	"""
	echo
	"""

}
