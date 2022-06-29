#!/usr/bin/env nextflow

include { download_bambi_decode_output_from_s3 } from '../ampseq-NF-pipeline/modules/download_test_data.nf'
include { bam_find_adapter } from '../ampseq-NF-pipeline/modules/bam_find_adapter.nf'
include { download_bamadapterfind_output_from_s3 } from '../ampseq-NF-pipeline/modules/download_test_data.nf'
include { compare_bam_subset } from '../ampseq-NF-pipeline/modules/test_tools.nf'

workflow {

	path_not_used = Channel.fromPath("./")
	test_ch = download_bambi_decode_output_from_s3("21045")
	input_bam = test_ch.map { it -> it[0] }

	test_input_ch = generate_test_channel(input_bam, path_not_used)

	under_test = bam_find_adapter(test_input_ch)
	test_bam = under_test.map { it -> it[1] }	

	reference_bam = download_bamadapterfind_output_from_s3("21045")

	compare_bam_subset(test_bam, reference_bam) 	
	 	



}



process generate_test_channel {

	input:
	path(bam_file)
	path(not_used1)

	output:
	tuple val("123_1"), path(bam_file), path(not_used1)

	script:
	"""
	echo
	"""

}
