#!/bin/bash

include { download_test_collated_bam_from_s3 } from  '../modules/download_test_data.nf'
include { alignment_filter } from '../modules/alignment_filter.nf'
include { check_md5sum } from '../modules/test_tools.nf'
include { check_md5sum as check_md5sum_metrics } from '../modules/test_tools.nf'

workflow {

	test_bam = download_test_collated_bam_from_s3("21045_1_ref")
	test_tag = "123"

	alignment_filter(test_tag, test_bam)		
	check_md5sum(alignment_filter.out.selected_bam, 'f70972064f62f3fa81750656e9c97391')
	check_md5sum_metrics(alignment_filter.out.selected_bam_metrics, 'e40015f9bbe0ab1d33cad9662081b53f')


}
