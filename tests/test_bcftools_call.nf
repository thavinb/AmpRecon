#!/usr/bin/env nextflow

include { download_ploidy_file_from_s3 } from '../modules/download_test_data.nf'
include { bcftools_call } from '../modules/bcftools_genotyping.nf'
include { download_test_file_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	reference_ploidy_file = download_ploidy_file_from_s3("Pf_GRC1v1.0")

	input_bcf_file = download_test_file_from_s3("21045_1_ref.mpileup.bcf")
	
	test_ch = create_tuple(input_bcf_file, reference_ploidy_file)
	bcftools_call(test_ch)
	called_file = bcftools_call.out.map{it -> it[1]}
	remove_header(called_file)
	check_md5sum(remove_header.out, "712687f4e782de9926f2e5d3c5263eec")
}

process remove_header {
    label 'bcftools'
    input:

    path(bcf_file)

    output:
    path("header_file")

    script:
    """
    bcftools view -H "${bcf_file}" > header_file
    """
}

process create_tuple {

    input:

    path(input_bcf)
    path(ploidy_file)

    output:
    tuple val("1"), path(input_bcf), path(ploidy_file)

    script:
    """
    echo 
    """


}
 

