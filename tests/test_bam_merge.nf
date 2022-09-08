#!/bin/bash

include { download_test_clipped_bam_from_s3; download_test_split_bam_from_s3 } from '../modules/download_test_data.nf'
include { bam_merge } from '../modules/bam_merge.nf'

workflow {

	//abandoned test - requires work re; container mount
	reheadered_bam = download_test_clipped_bam_from_s3("21045_1_1")
	split_bam = download_test_split_bam_from_s3("21045_1_1")
	test_tag = "123"

	test_ch = generate_test_channel(reheadered_bam, split_bam, test_tag)
	bam_merge(test_ch)


}


process generate_test_channel {

	input:
	path(reheadered_bam)
	path(split_bam)
	val(test_tag)

	output:
	tuple val(test_tag), path(reheadered_bam), path(split_bam)

	script:
	"""
	echo
	"""

}
