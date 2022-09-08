#!/usr/bin/env nextflow

include { download_bambi_decode_output_from_s3 } from '../modules/download_test_data.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { download_bamadapterfind_output_from_s3 } from '../modules/download_test_data.nf'
include { download_i2b_output_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	path_not_used = Channel.fromPath("./")
	test_ch = download_i2b_output_from_s3("21045")

	test_input_ch = generate_test_channel(test_ch, path_not_used)

	under_test = bam_find_adapter(test_input_ch)
	test_bam = under_test.map { it -> it[1] }	


	check_md5sum(test_bam, "5153ba386b93425fc1066ef7c66e9b6f") 	
	 	



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
