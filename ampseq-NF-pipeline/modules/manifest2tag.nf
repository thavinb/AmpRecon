#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process get_tag_list_file {
	input:
	path manifest
	val library
	val sample
	val study

	output:
	path("${tag_list}"), emit: taglist_file

	script:
	tag_list = "tag_file.tsv"
	"""
	python ${workflow.projectDir}/modules/manifest2tag.py -m $manifest -l $library --sample $sample --study $study
	"""

}
