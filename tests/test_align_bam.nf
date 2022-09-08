#!/usr/bin/env nextflow

include { download_test_fastq_from_s3 } from '../modules/download_test_data.nf'
include { download_idx_from_s3 } from '../modules/download_test_data.nf'
include { download_fasta_from_s3 } from '../modules/download_test_data.nf'
include { align_bam } from '../modules/align_bam.nf'
include { download_test_aligned_sam_from_s3 } from '../modules/download_test_data.nf'
include { test_size } from '../modules/test_tools.nf'
include { test_binary } from '../modules/test_tools.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	index = download_idx_from_s3("Pf_GRC1v1.0")
	ref = download_fasta_from_s3("Pf_GRC1v1.0")
	fasta = ref.map { it -> it[0] }
	fastq = download_test_fastq_from_s3("21045_1_ref")
	sample_tag = "1"
	run_id = "21045"

	align_bam(sample_tag, fastq, fasta, index, run_id)
	
	test_sam= align_bam.out.sam_file

	check_md5sum(test_sam, "6b85be31d45052fd57ec8f4147b045e9") 
		

	



}

