#!/usr/bin/env nextflow

include { download_bambi_decode_output_from_s3 } from '../modules/download_test_data.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { download_bamadapterfind_output_from_s3 } from '../modules/download_test_data.nf'
include { download_test_cram_from_s3 } from '../modules/download_test_data.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'
include { check_cram_md5sum } from '../modules/test_tools.nf'


workflow {

	run_id = "21045"
	download_bamadapterfind_output_from_s3("21045")
	test_bam = download_bamadapterfind_output_from_s3.out.test_bam
	test_metrics = download_bamadapterfind_output_from_s3.out.metrics	

	under_test = bam_to_cram(run_id, test_bam, test_metrics)
	cram = under_test.cram_fls
	
	check_cram_md5sum(cram, "06ffc3ed4974a1a34c39485bbef1381f") 
		
}


