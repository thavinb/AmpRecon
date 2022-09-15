#!/usr/bin/env nextflow


include { download_i2b_output_from_s3 } from '../modules/download_test_data.nf'
include { download_bambi_decode_output_from_s3 } from '../modules/download_test_data.nf'
include { decode_multiplexed_bam } from '../modules/decode_multiplexed_bam.nf'
include { check_md5sum } from '../modules/test_tools.nf'
include { check_md5sum as check_md5sum_metrics } from '../modules/test_tools.nf'

workflow {
	
	input_bam = download_i2b_output_from_s3("21045")
	input_tags = Channel.fromPath("./tags.tsv")
	base_dir = Channel.fromPath("./")	

	decode_ch = generate_test_channel(input_bam, input_tags, base_dir)

	under_test = decode_multiplexed_bam(decode_ch)
	test_bam = under_test.map { it -> it[1] }
	test_metrics = under_test.map { it -> it[2] }

	reference_data = download_bambi_decode_output_from_s3("21045")
	reference_bam = reference_data.map { it -> it[0] }
	reference_metrics = reference_data.map { it -> it[1] }

	test_preprocessing = remove_id_tags_in_metrics_header(reference_metrics, test_metrics)
	reference_metrics = test_preprocessing.map { it -> it[0] }
	test_metrics = test_preprocessing.map { it -> it[1] }

	check_md5sum(test_bam, "e3f44394583b51fdefd0d8e90c0a2d09")
	check_md5sum_metrics(test_metrics, "0852211f9d8e476dabf4772a3ba51306")
	 


}


process generate_test_channel {
        //to circumvent a lot of the pain around getting correct input channels for bambi i2b process this stub process outputs
        //the channel in correct fashion

        input:
        path(input_bam)
        path(tags_file)
	path(base_dir) //for the barcodes as input

        output:
        tuple val("21045"), path(base_dir),  val("1"), val("Test"), val("21045"), val("lib"), path(tags_file), path(input_bam)

        script:
        """
        echo
        """
}

process remove_id_tags_in_metrics_header {


	input:
	path(reference_metrics)
	path(test_metrics)

	output:
	tuple path(reference_output), path(test_output)

	script:
	reference_output="reference.metrics"
	test_output="test.metrics"
	"""
	grep -v "ID:" ${reference_metrics} > ${reference_output}
	grep -v "ID:" ${test_metrics} > ${test_output}
	"""

}

