#!/bin/bash

include { download_test_clipped_bam_from_s3; download_test_scrambled_bam_from_s3;
	  download_fasta_from_s3	} from '../modules/download_test_data.nf'
include { mapping_reheader } from '../modules/mapping_reheader.nf' 


workflow {

	test_clipped_bam = download_test_clipped_bam_from_s3("21045_1_1")
	test_scrambled_bam = download_test_scrambled_bam_from_s3("21045_1_1")
	test_tag = "123"

	test_ch = generate_test_channel(test_clipped_bam, test_scrambled_bam, test_tag)

	reference = download_fasta_from_s3("21045_1_1")
	reference_fasta = reference.map { it -> it[0] }
	reference_dict = reference.map { it -> it[1] }

	mapping_reheader(test_ch, reference_fasta, reference_dict)


}


process generate_test_channel {


	input:
	path(clipped_bam)
	path(scrambled_bam)
	val(test_tag)

	output:
	tuple val(test_tag), path(scrambled_bam), path(clipped_bam)

	script:
	"""
	echo
	"""
}

