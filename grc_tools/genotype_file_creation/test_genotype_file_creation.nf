#!/usr/bin/env nextflow

include { assemble_genotype_file } from '../modules/assemble_genotype_file.nf'
include { download_test_file_from_s3 } from '../modules/download_test_data.nf'
include { download_test_file_from_s3 as download_first_vcf } from '../modules/download_test_data.nf'
include { download_test_file_from_s3 as download_second_vcf } from '../modules/download_test_data.nf'
include { download_test_file_from_s3 as download_third_vcf } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	input_vcf_file_1 = download_first_vcf("21045_1_ref.bcftools_genotyped.vcf.gz")
        input_vcf_file_2 = download_second_vcf("21045_2_ref.bcftools_genotyped.vcf.gz")
        chromKey_file = download_test_file_from_s3("chromKey.txt")

	//test_ch = create_tuple(input_vcf_file_1, input_vcf_file_2)
        test_ch = input_vcf_file_1.concat(input_vcf_file_2).map{it -> tuple("21045", it)}.groupTuple()
        assemble_genotype_file(test_ch, chromKey_file)
        genotype_file = assemble_genotype_file.out.map{it -> it[1]}

        // Check md5 hash of the genotypes file
	check_md5sum(genotype_file, "469b9753db48d71b539ff9f3eff2cd8a")
}

process create_tuple {

    input:

    val(input_vcf_1)
    val(input_vcf_2)

    output:
    tuple val("1"), val(file_list)

    script:
    file_list = "${input_vcf_1} ${input_vcf_2}"
    """
    echo 
    """
}

