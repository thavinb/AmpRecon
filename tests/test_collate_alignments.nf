#!/usr/bin/env nextflow

include { download_test_cram_from_s3 } from '../modules/download_test_data.nf'
include { download_test_collated_bam_from_s3 } from '../modules/download_test_data.nf'
include { collate_alignments } from '../modules/collate_alignments.nf'
include { check_md5sum } from '../modules/test_tools.nf'


workflow {

	input_cram = download_test_cram_from_s3("21045_1_ref")
	test_ch = generate_test_channel(input_cram)

	under_test = collate_alignments(test_ch)
	test_bam = under_test.collated_bam

	check_md5sum(test_bam, "f290d2d1c1a1794b4e45e611e1dd4c57")
}


process generate_test_channel {


	input:
	path(input_cram)

	output:
	val(21045)
        path(input_cram)
        val("123_1")

	script:
	"""
	echo
	"""
}
