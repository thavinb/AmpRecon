#!/usr/bin/env nextflow

include { download_annotation_vcf_from_s3 } from '../modules/download_test_data.nf'
include { download_idx_from_s3 } from '../modules/download_test_data.nf'
include { download_fasta_from_s3 } from '../modules/download_test_data.nf'
include { bcftools_mpileup } from '../modules/bcftools_genotyping.nf'
include { download_test_file_from_s3; download_test_file_from_s3 as download_another_test_file_from_s3 } from '../modules/download_test_data.nf'
include { check_md5sum } from '../modules/test_tools.nf'

workflow {

	reference_index = download_idx_from_s3("Pf_GRC1v1.0")
	download_fasta_from_s3("Pf_GRC1v1.0").map{ it -> it[0] }.set{reference}
	reference_annotation = download_annotation_vcf_from_s3("Pf_GRC1v1.0")

	download_test_file_from_s3("21045_1_ref.samtools_sorted.bam")
	download_another_test_file_from_s3("21045_1_ref.samtools_sorted.bam.bai")

	bam = download_test_file_from_s3.out
	bam_index = download_another_test_file_from_s3.out

	test_ch = create_tuple(bam, bam_index, reference, reference_index, reference_annotation)
	bcftools_mpileup(test_ch)
	mpileup_file = bcftools_mpileup.out.map{it -> it[1]}
	check_md5sum(mpileup_file, "c49b45f99569631a33abaf49899cd35c")
}


process create_tuple {

    input:

    path(bam)
    path(bam_index)
    path(ref)
    path(ref_index)
    path(ref_annotation)

    output:
    tuple val("1"), path(bam), path(bam_index), path(ref), path(ref_index), path(ref_annotation)

    script:
    """
    echo 
    """


}
 

