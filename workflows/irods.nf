#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import subworkflows
include { PULL_FROM_IRODS } from './pipeline-subworkflows/pull_from_irods.nf'


workflow IRODS {
    take:
        irods_manifest // irods manifest file
        reference_ch // tuple ([fasta], panel_name, [fasta_idx_files], [dictionary_file], [ploidy_file], [annotation_vcf_file], [snp_list])
    main:
        // load manifest content
        irods_ch =  Channel.fromPath(irods_manifest, checkIfExists: true)
                      | splitCsv(header: true, sep: '\t')
                      //| map { row -> tuple(row.id_run, row.primer_panel, row.WG_lane) }
                      | map { row ->
                              WG_lane = "${row.irods_path}".split('/')[-1].split('\\.')[0]
                              tuple(row.sample_id, row.primer_panel, WG_lane, row.irods_path) 
                            }
//reference_ch.first().view()
        // Assign each sample id the appropriate set of reference files
        irods_ch
             | combine(reference_ch,  by: 1) // tuple (primer_panel, sample_id, WG_lane, irods_path, [fasta], [fasta_idx_files], [dictionary_file], [ploidy_file], [annotation_vcf_file], [snp_list])
             | map { it -> tuple(it[2], it[1], it[4], it[5], it[0], it[6], it[7], it[8], it[9]) }
             | set{ sample_id_ref_ch } // tuple (WG_lane, sample_id, fasta_file, fasta_idx, primer_panel, [dictionary_file], [ploidy_file], [annotation_vcf_file], [snp_list])
// sample_id_ref_ch.first().view() //--checkpoint
        // remove panels info from channel (is not used on this subworkflow)
        irods_ch.map{ it -> tuple (it[0], it[2], it[3]) }.set{irods_ch_noRef} // tuple(sample_id, WG_lane, irods_path)

        // run step1.2b - pull from iRODS
        PULL_FROM_IRODS(irods_ch_noRef, sample_id_ref_ch)

        // prepare channel for step 1.3
        sample_tag_reference_files_ch = PULL_FROM_IRODS.out.sample_tag_reference_files_ch

        bam_files_ch = PULL_FROM_IRODS.out.bam_files_ch
    emit:
        bam_files_ch //irods_Out_ch
        sample_tag_reference_files_ch
}



