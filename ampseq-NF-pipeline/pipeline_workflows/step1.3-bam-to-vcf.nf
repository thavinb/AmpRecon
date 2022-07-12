#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

include { bam_reset } from '../modules/bam_reset.nf'


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
