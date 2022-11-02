#!/usr/bin/env nextflow

include { bcftools_filter } from '../modules/bcftools_genotyping.nf'
include { download_test_file_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	input_bcf_file = download_test_file_from_s3("21045_1_ref.call.bcf")
	
	test_ch = create_tuple(input_bcf_file)
	bcftools_filter(test_ch)
	filtered_file = bcftools_filter.out.map{it -> it[1]}
        remove_header(filtered_file)
	check_md5sum(remove_header.out, "8c6d7566cc55726baccdba916be28d9d")
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

    output:
    tuple val("1"), path(input_bcf)

    script:
    """
    echo 
    """


}
 

