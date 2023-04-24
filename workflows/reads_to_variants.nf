#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { ALIGNMENT } from './alignment.nf'
include { READ_COUNTS } from './read_counts.nf'
include { GENOTYPING_GATK } from './genotyping_gatk.nf'
include { GENOTYPING_BCFTOOLS } from './genotyping_bcftools.nf'
include { write_vcfs_manifest } from '../modules/write_vcfs_manifest.nf'
/*
Here all workflows which are used regardless of the entry point (iRODS or inCountry)
are setup
*/

workflow READS_TO_VARIANTS {
    take:
        fastq_ch // tuple (file_id, fastq_file, reference_fasta_file)
        file_id_reference_files_ch // tuple (file_id, panel_name, reference_fasta_file, snp_list)
        annotations_ch // tuple (panel_name, annotation_file)
        file_id_to_sample_id_ch // tuple (file_id, sample_id)

    main:

        // alignment
        file_id_reference_files_ch.map{it -> tuple(it[0], it[2], it[1])}.set{alignment_ref_ch} // tuple (file_id, fasta_file, panel_name)
        ALIGNMENT(fastq_ch)
        
	//read counts 
	READ_COUNTS(ALIGNMENT.out, alignment_ref_ch, annotations_ch)

        // genotyping
        if( params.genotyping_gatk == true ) {
            file_id_reference_files_ch.map{it -> tuple(it[0], it[2], it[3])}.set{gatk_genotyping_ref_ch} // tuple (file_id, fasta_file, snp_list)
            GENOTYPING_GATK(
                ALIGNMENT.out,
                gatk_genotyping_ref_ch
            )
            GENOTYPING_GATK.out.set{genotyping_gatk_ch}
        } else {
            Channel.empty().set{ genotyping_gatk_ch }
        }

        if( params.genotyping_bcftools == true ) {
            file_id_reference_files_ch.map{it -> tuple(it[0], it[2], it[3])}.set{bcftools_genotyping_ref_ch} // tuple (file_id, fasta_file, snp_list)
            GENOTYPING_BCFTOOLS(
                ALIGNMENT.out,
                bcftools_genotyping_ref_ch
            )
            GENOTYPING_BCFTOOLS.out.set{genotyping_bcftools_ch}
        } else {
            Channel.empty().set{ genotyping_bcftools_ch }
        }

        // concatonate genotyping workflow output channels
        genotyping_bcftools_ch.concat(genotyping_gatk_ch).set{vcf_files_ch}

        // Create channel of 2 lists: IDs and VCFs
        vcf_files_ch // tuple (file_id, vcf_path, vcf_index_path)
            .map{it -> tuple(it[0], it[1])}
            .join(file_id_to_sample_id_ch) // tuple (file_id, vcf_path, sample_id)
            .multiMap { it ->
                id_list: it[2]
                vcf_list: it[1]
            }.set{lanelet_ch}
    
        // Write manifest of lanelet VCF files and sample IDs
        write_vcfs_manifest(lanelet_ch.id_list.collect(), lanelet_ch.vcf_list.collect())
        lanelet_manifest = write_vcfs_manifest.out

    emit:
        lanelet_manifest
}

