#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { grc_kelch13_mutation_caller } from '../grc_tools/kelch13/grc_kelch13_mutation_caller.nf'
include { grc_plasmepsin_cnv_caller } from '../grc_tools/plasmepsin/grc_plasmepsin_cnv_caller.nf'
include { grc_speciate } from '../grc_tools/speciation/grc_speciate.nf'
include { grc_barcoding } from '../grc_tools/barcode/grc_barcoding.nf'
include { grc_estimate_coi } from '../grc_tools/COI/grc_estimate_coi.nf'
include { grc_amino_acid_caller } from '../grc_tools/amino_acid_calling/grc_amino_acid_caller.nf'
include { grc_assemble } from '../grc_tools/assemble_grc1/grc_assemble.nf'

workflow GENOTYPES_TO_GRCS {
    take:
        genotype_files_ch
        kelch_reference_file
        codon_key_file
        drl_information_file

    main:
        // Call mutations at Kelch13 loci
        grc_kelch13_mutation_caller(genotype_files_ch, kelch_reference_file, codon_key_file)

        // Call copy number variation at Plasmepsin breakpoint
        grc_plasmepsin_cnv_caller(genotype_files_ch)
        
        // Create barcodes
        grc_barcoding(genotype_files_ch)

        // Determine species
        grc_speciate(genotype_files_ch, grc_barcoding.out)

        // Complexity of infection estimation
        grc_estimate_coi(grc_barcoding.out)

        // Assemble drug resistance haplotypes and GRC2
        grc_amino_acid_caller(genotype_files_ch, drl_information_file, codon_key_file)

        // Assemble GRC1
        grc_kelch13_mutation_caller.out
            .concat(grc_plasmepsin_cnv_caller.out)
            .concat(grc_barcoding.out)
            .concat(grc_speciate.out)
            .concat(grc_estimate_coi.out)
            .concat(grc_amino_acid_caller.out.drl_haplotypes)
            .collect()
            .set{grc1_columns}
        grc_assemble(grc1_columns)
}

workflow {
    genotype_files_ch = Channel.fromPath(params.genotype_files_path, checkIfExists: true).collect()
    kelch_reference_file = Channel.fromPath(params.kelch_reference_file_path, checkIfExists: true)
    codon_key_file = Channel.fromPath(params.codon_key_file_path, checkIfExists: true)
    drl_information_file = Channel.fromPath(params.drl_information_file, checkIfExists: true)

    GENOTYPES_TO_GRCS(genotype_files_ch, kelch_reference_file, codon_key_file, drl_information_file)
}

