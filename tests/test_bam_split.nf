#!/bin/bash

include { download_test_reheadered_bam_from_s3 } from '../modules/download_test_data.nf'
include { bam_split } from '../modules/bam_split.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	test_tag = "123"
	reheadered_bam = download_test_reheadered_bam_from_s3("21045_1_ref")

	test_ch = generate_test_ch(test_tag, reheadered_bam)

	under_test = bam_split(test_ch)
	bam = under_test.map { it -> it[1] }

	check_md5sum(bam, "23915452b9c7dfeca6ae411b9deb87b2") 


}


process generate_test_ch {

	input:
	val(test_tag)
	path(reheadered_bam)


	output:
	tuple val(test_tag), path(reheadered_bam)

	script:
	"""
	echo 
	"""
}
