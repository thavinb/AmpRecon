#!/usr/bin/env nextflow

include { subset_bam } from '../ampseq-NF-pipeline/modules/test_tools.nf'
include { compare_bam_subset } from '../ampseq-NF-pipeline/modules/test_tools.nf'
include { test_binary } from '../ampseq-NF-pipeline/modules/test_tools.nf'
include { download_redo_alignment_output_bam_from_s3 } from '../ampseq-NF-pipeline/modules/download_test_data.nf'
include { download_zipped_read_counts_from_s3 } from '../ampseq-NF-pipeline/modules/download_test_data.nf'

nextflow.enable.dsl = 2

workflow { 
        /*
        Validation of the output files produced by a test run 
        of the pipeline in-country entry point.
        */
        
        // Validation of Step1.3 BAM files

	test_bam = Channel.fromPath("${params.results_dir}/27122/27122_1#100.bam")
        reference_bam = download_redo_alignment_output_bam_from_s3("27122_1#100")

	//test_bam_subset = subset_bam(test_bam) //subset to make more efficient test for differences
	compare_bam_subset(test_bam, reference_bam)//test_bam_subset, reference_bam) 

        // Validation of Step1.3 read count files
        read_count_files = Channel.fromPath( [file("${params.results_dir}/27122/27122_GRC1_reads_per_region.csv"), 
                                              file("${params.results_dir}/27122/27122_GRC2_reads_per_region.csv"), 
                                              file("${params.results_dir}/27122/27122_Spec_reads_per_region.csv")] )

        download_zipped_read_counts_from_s3("27122")
        unzip_read_counts(download_zipped_read_counts_from_s3.out)
        reference_read_count_files = unzip_read_counts.out.flatten()
        
        test_binary(reference_read_count_files, read_count_files)
}

process unzip_read_counts {


	input:
            path(input_read_counts_zip)

	output:
            file("*_reads_per_region.csv")
            //path("*_reads_per_region.csv")
	script:

	"""
        tar xvzf "${input_read_counts_zip}" -C .
	"""

}

