#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { REALIGNMENT } from './pipeline-subworkflows/realignment.nf'
include { GENOTYPING } from './pipeline-subworkflows/genotyping.nf'

/*
Here all workflows which are used regardless of the entry point (iRODS or inCountry)
are setup
*/

workflow COMMON {

    take:
        bam_files_ch // tuple(sample_id, bam_file)
        sample_tag_reference_files_ch // tuple(sample_id, ref_files)
        annotations_ch // tuple (pannel_name, anotation_file)
    main:
        // mapping tuple to multichannel 
        bam_files_ch
            | multiMap {
                sample_tag: it[0]
                bam_file: it[1]
                }
            | set { realignment_In_ch }
        
        // do realignment and read counts
        REALIGNMENT(
                    realignment_In_ch.sample_tag,
                    realignment_In_ch.bam_file,
                    sample_tag_reference_files_ch,
                    annotations_ch
                )
        
        // genotyping
        GENOTYPING(
                   REALIGNMENT.out,
                   sample_tag_reference_files_ch
                )
}
