#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { REALIGNMENT } from './pipeline-subworkflows/realignment.nf'
include { GENOTYPING_GATK } from './pipeline-subworkflows/genotyping_gatk.nf'
include { GENOTYPING_BCFTOOLS } from './pipeline-subworkflows/genotyping_bcftools.nf'

include { bqsr } from '../modules/bqsr.nf' addParams(gatk:params.gatk3)
include { samtools_index } from '../modules/samtools.nf'
/*
Here all workflows which are used regardless of the entry point (iRODS or inCountry)
are setup
*/

workflow COMMON {

    take:
        bam_files_ch // tuple(sample_id, bam_file)
        sample_tag_reference_files_ch // tuple(sample_id, ref_files)
        annotations_ch // tuple (panel_name, anotation_file)
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
        
        BQSR(
            REALIGNMENT.out,
            sample_tag_reference_files_ch
        )

        // genotyping
        if( params.genotyping_gatk == true ) {
                GENOTYPING_GATK(
                   BQSR.out,
                   sample_tag_reference_files_ch
                )
        }

        if( params.genotyping_bcftools == true ) {
                GENOTYPING_BCFTOOLS(   
                   BQSR.out,
                   sample_tag_reference_files_ch
                )
        }
}

workflow BQSR {
    take:
        input_sample_tags_bams_indexes
        sample_tag_reference_files_ch
    
    main:
        // base quality score recalibration
        input_sample_tags_bams_indexes // tuple (sample_tag, bam_file, bam_index)
            | join(sample_tag_reference_files_ch) // tuple (sample_tag, bam_file, bam_index [ref_files])
            | map{ it -> tuple(it[0], it[1], it[2], it[3], it[3]+".fai", it[3]+".dict")}
            | set{bqsr_input} // tuple(sample_tag, bam_file, bam_index, [ref_files], )

        if (!params.skip_bqsr) {
            bqsr(bqsr_input)
            samtools_index(bqsr.out)
            
            // haplotype caller
            post_bqsr_output = samtools_index.out
        }
        else {
            post_bqsr_output = input_sample_tags_bams_indexes
        }
    emit:
        post_bqsr_output
}