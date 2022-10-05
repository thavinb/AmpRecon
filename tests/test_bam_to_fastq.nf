#!/usr/bin/env nextflow

include { download_test_clipped_bam_from_s3 } from '../modules/download_test_data.nf'
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'
include { download_test_fastq_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	test_bam = download_test_clipped_bam_from_s3("21045_1_ref")
	test_sample_tag="123_1"

	under_test = bam_to_fastq(test_sample_tag, test_bam)
	test_fastq = under_test.map { it -> it[1] }

	check_md5sum(test_fastq, "869cb8be189b0a4b57f2ba634ff575a6") 

}
