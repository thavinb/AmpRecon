#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { ALIGNMENT } from './alignment.nf'
include { READ_COUNTS } from './read_counts.nf'
include { GENOTYPING_GATK } from './genotyping_gatk.nf'
include { GENOTYPING_BCFTOOLS } from './genotyping_bcftools.nf'
include { write_vcfs_manifest } from '../modules/write_vcfs_manifest.nf'
include { merge_bams_and_index } from '../modules/bam_merge_and_index.nf'
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
			   
	ALIGNMENT.out
	         .join( file_id_to_sample_id_ch )
		 .join ( file_id_reference_files_ch )
		 .map { it -> tuple("${it[3]}_${it[4]}", it[1]) } 
		 .groupTuple()
		 .set { bam_merge_ch } 
 
        // read counts
        READ_COUNTS(ALIGNMENT.out, alignment_ref_ch, annotations_ch)

        merge_bams_and_index( bam_merge_ch )

	file_id_reference_files_ch.join( file_id_to_sample_id_ch )
				  .map { it -> tuple("${it[4]}_${it[1]}", it[0], it[1], it[2], it[3], it[4] ) }
				  .set {test_ch}
	

        // genotyping
        if( params.genotyping_gatk == true ) {
            test_ch.map{it -> tuple(it[0], it[1], it[3], it[4])}.set{gatk_genotyping_ref_ch} // tuple (file_id, fasta_file, snp_list, sample_key)
            GENOTYPING_GATK(
                merge_bams_and_index.out,
                gatk_genotyping_ref_ch
            )
            GENOTYPING_GATK.out.set{genotyping_gatk_ch}
        } else {
            Channel.empty().set{ genotyping_gatk_ch }
        }

        if( params.genotyping_bcftools == true ) {
            test_ch.map{it -> tuple(it[0],it[1], it[3], it[4])}.set{bcftools_genotyping_ref_ch} // tuple (file_id, fasta_file, snp_list, sample_key)
            GENOTYPING_BCFTOOLS(
                merge_bams_and_index.out,
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

