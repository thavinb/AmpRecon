#!/bin/bash

include { download_test_collated_bam_from_s3 } from '../modules/download_test_data.nf'
include { sort_bam } from '../modules/sort_bam.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	test_bam = download_test_collated_bam_from_s3("21045_1_ref")
	test_run_id = "2001"
	test_tag = "123"

	under_test = sort_bam(test_run_id, test_tag, test_bam)
	
	output_bam = under_test.map { it -> it[1] }

	check_md5sum(output_bam, "56589d8ecc063fef41ac1b0dd4e06287")



}
