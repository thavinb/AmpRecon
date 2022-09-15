#!/usr/bin/env nextflow

include { download_bambi_decode_output_from_s3 } from '../modules/download_test_data.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { download_bamadapterfind_output_from_s3 } from '../modules/download_test_data.nf'
include { download_test_cram_from_s3 } from '../modules/download_test_data.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'
include { compare_bam_subset } from '../modules/test_tools.nf'


workflow {

	input_bam = download_bamadapterfind_output_from_s3("21045")
	test_ch = generate_test_channel(input_bam)

	under_test = bam_to_cram(test_ch)
	test_bam = under_test.map { it -> it[1] }
	reference_bam = download_test_cram_from_s3("21045_1_ref")
	
	compare_bam_subset(test_bam, reference_bam)
}


process generate_test_channel {


	input:
	path(input_bam)

	output:
	tuple val(21045), path(input_bam)

	script:
	"""
	echo
	"""
}

