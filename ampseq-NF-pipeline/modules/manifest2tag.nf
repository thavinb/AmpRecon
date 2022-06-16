#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process get_taglist_file {
	publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true
	input:
	tuple val(run_id), path(bcl_dir), val(lane), val(study_name), val(read_group), val(library), val(manifest)
	//path manifest
	//val library
	//val sample
	//val study

	output:
	path("${tag_list}"), emit: taglist_file

	script:
	tag_list = "tag_file.tsv"
	// not sure if run id shoud go here, check 210714.taglist as reference
	// current understanding is that this information is not used
  sample = run_id
	"""
	python ${workflow.projectDir}/modules/manifest2tag.py -m $manifest -l $library --sample $sample --study $study_name
	"""

}
