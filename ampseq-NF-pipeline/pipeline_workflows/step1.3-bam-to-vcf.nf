#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../modules/bam_reset.nf'
include { align_bam } from '../modules/align_bam.nf'
include {indexReferenceBWA} from "../modules/bwa_index.nf"
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'


workflow reset_bam_alignment {

  //remove alignment from bam - this process proceeds directly after the end of 1.2

  take:
    //input_ch // [ sample_tag, bam_file, run_id ]
    sample_tag
    bam_file
    run_id
    //input_manifest //from step 1.2, path/file object what about irods

  main:
    /*input_ch = input_manifest.splitCsv(header : true)
              .multiMap {
                         row  -> run_id:row.run_id
                                 bam_fl:row.bam_fl
                                 sample_tag:row.sample_tag
              }
    */
    bam_reset(sample_tag, bam_file)
}

workflow read_alignment {

	take:
		reference_genome
		input_fastq_ch //this can be output from bam_to_fastq


	main:
		sample_tag = input_fastq_ch.map { it -> it[0] } 
		input_fastq = input_fastq_ch.map { it -> it[1] }
		reference_idx = indexReferenceBWA(reference_genome) //could add this as input to the workflow
		align_bam(sample_tag, input_fastq, reference_genome, reference_idx) 

}

