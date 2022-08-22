include { pull_from_iRODS } from './pipeline-subworkflows/step1.2b-pull-from-iRODS.nf'
include { redo_alignment } from './pipeline-subworkflows/step1.3-redo_alignment.nf'

workflow IRODS {
	
	take:
	input_params_file // manifest/input.csv

	main:
	pull_from_iRODS(input_params_file) //get crams
	redo_alignment(pull_from_irods.out)

	//stage two here

}



