#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { REALIGNMENT } from './pipeline-subworkflows/realignment.nf'
include { GENOTYPING_GATK } from './pipeline-subworkflows/genotyping_gatk.nf'
include { GENOTYPING_BCFTOOLS } from './pipeline-subworkflows/genotyping_bcftools.nf'

/*
Here all workflows which are used regardless of the entry point (iRODS or inCountry)
are setup
*/

workflow COMMON {

    take:
        bam_files_ch // tuple(sample_id, bam_file)
        sample_tag_reference_files_ch // tuple(new_sample_id, fasta_file, ref_files, fasta_index, panel_names, [dictionary_file], [ploidy_file], [annotation_vcf_file], [snp_list])
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
        sample_tag_reference_files_ch.map{it -> tuple(it[0], it[1], it[2], it[3])}.set{realignment_input_reference_ch}
        REALIGNMENT(
                    realignment_In_ch.sample_tag,
                    realignment_In_ch.bam_file,
                    realignment_input_reference_ch,
                    annotations_ch
                )
        realignment_input_reference_ch.first().view()
        // genotyping
        if( params.genotyping_gatk == true ) {
        sample_tag_reference_files_ch.map{it -> tuple()}.set{gatk_genotyping_input_reference_ch}
                GENOTYPING_GATK(
                   REALIGNMENT.out,
                   gatk_genotyping_input_reference_ch
                )
        }

        if( params.genotyping_bcftools == true ) {
        sample_tag_reference_files_ch.map{it -> tuple()}.set{bcftools_genotyping_input_reference_ch}
                GENOTYPING_BCFTOOLS(   
                   REALIGNMENT.out,
                   bcftools_genotyping_input_reference_ch
                )
        }
}
