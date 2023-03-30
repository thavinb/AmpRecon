#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { ALIGNMENT } from './alignment.nf'
include { GENOTYPING_GATK } from './genotyping_gatk.nf'
include { GENOTYPING_BCFTOOLS } from './genotyping_bcftools.nf'
include { samtools_index } from '../modules/samtools.nf'
include { write_vcfs_manifest } from '../modules/write_vcfs_manifest.nf'
/*
Here all workflows which are used regardless of the entry point (iRODS or inCountry)
are setup
*/

workflow READS_TO_VARIANTS {

    take:
        bam_files_ch // tuple(sample_tag, bam_file)
        sample_tag_reference_files_ch // tuple (sample_tag, panel_name, reference_fasta_file, snp_list)
        annotations_ch // tuple (panel_name, annotation_file)
        file_id_to_sample_id_ch // tuple (sample_tag, sample_id)
    main:
        // mapping tuple to multichannel 
        if (params.aligned_bams_mnf == null){
            bam_files_ch
              | multiMap {
                sample_tag: it[0]
                bam_file: it[1]
                }
              | set { alignment_In_ch }
        
        // do alignment and read counts
        sample_tag_reference_files_ch.map{it -> tuple(it[0], it[2], it[1])}.set{alignment_ref_ch} // tuple (sample_tag, fasta_file, panel_name)
        ALIGNMENT(
                    alignment_In_ch.sample_tag,
                    alignment_In_ch.bam_file,
                    alignment_ref_ch,
                    annotations_ch
                )

        genotyping_In_ch = ALIGNMENT.out
        }

        if (params.execution_mode == "aligned_bams"){
        genotyping_In_ch = bam_files_ch
        }
        // genotyping
        sample_tag_reference_files_ch.map{it -> tuple(it[0], it[2], it[3])}.set{bqsr_ref_ch} // tuple (sample_tag, fasta_file, snp_list)
        BQSR(
            genotyping_In_ch,
            bqsr_ref_ch
        )

        // genotyping
        if( params.genotyping_gatk == true ) {
        sample_tag_reference_files_ch.map{it -> tuple(it[0], it[2], it[3])}.set{gatk_genotyping_ref_ch} // tuple (sample_tag, fasta_file, snp_list)
                GENOTYPING_GATK(
                   BQSR.out,
                   gatk_genotyping_ref_ch
                )
                GENOTYPING_GATK.out.set{genotyping_gatk_ch}
        } else {
            Channel.empty().set{ genotyping_gatk_ch }
        }

        if( params.genotyping_bcftools == true ) {
        sample_tag_reference_files_ch.map{it -> tuple(it[0], it[2], it[3])}.set{bcftools_genotyping_ref_ch} // tuple (sample_tag, fasta_file, snp_list)
                GENOTYPING_BCFTOOLS(
                   BQSR.out,
                   bcftools_genotyping_ref_ch
                )
                GENOTYPING_BCFTOOLS.out.set{genotyping_bcftools_ch}
        } else {
            Channel.empty().set{ genotyping_bcftools_ch }
        }

        // concatonate genotyping workflow output channels
        genotyping_bcftools_ch.concat(genotyping_gatk_ch).set{vcf_files_ch}

        // Prepare manifest of lanelet VCF files and sample IDs for GRC creation
        vcf_files_ch // tuple (file_id, vcf_path, vcf_index_path)
            .map{it -> tuple(it[0], it[1])}
            .join(file_id_to_sample_id_ch) // tuple (file_id, vcf_path, sample_id)
            .multiMap { it ->
                id_list: it[2]
                vcf_list: it[1] }
            .set{lanelet_ch}
        write_vcfs_manifest(lanelet_ch.id_list.collect(), lanelet_ch.vcf_list.collect())
        lanelet_manifest = write_vcfs_manifest.out

    emit:
        lanelet_manifest
}

