#!/usr/bin/env nextflow

include { subset_bam } from '../modules/test_tools.nf'
include { compare_bam_subset } from '../modules/test_tools.nf'
include { basecalls_conversion } from '../modules/basecalls_conversion.nf'
include { download_i2b_output_from_s3 } from '../modules/download_test_data.nf'
include { download_bcl_from_s3 } from '../modules/download_test_data.nf'

nextflow.enable.dsl = 2


workflow {

	// container should be included later

	bcl_dir_ch = download_bcl_from_s3("21045") //gets test bcl directory from s3
	tags_file = Channel.fromPath("./tags.tsv") //fictitious data - not sure what purpose it serves

	i2b_input_ch = generate_test_channel(bcl_dir_ch, tags_file) //under test 
	test_bam = basecalls_conversion(i2b_input_ch)
	
	reference = download_i2b_output_from_s3("21045")
	under_test = subset_bam(test_bam) //subset to make more efficient test for differences
	
	compare_bam_subset(under_test, reference) 

}


process generate_test_channel {
	//to circumvent a lot of the pain around getting correct input channels for bambi i2b process this stub process outputs
	//the channel in correct fashion	

	input:
	path(bcl_dir)
	path(tags_file)

	output:
	tuple val("21045"), path(bcl_dir), val("1"), val("Test"), val("21045"), val("lib"), path(tags_file)
	
	script:
	"""
	echo  
	"""
}


