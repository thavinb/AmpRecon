#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {indexReferenceBWA} from "../modules/bwa_index.nf"
include {generateFastaIndex} from "../modules/generateFastaIndex.nf"
include {generateFastaDict} from "../modules/generateFastaDict.nf"

workflow prepare_reference {
  take:
    reference_fasta // <string> the path to a reference fasta file

  main:
    // --- Sanity Check -------------------------------------------------------
    // check if fasta exists and follow symlinks if needed
    //ref_fa = file(reference_fasta, checkIfExists=true, followLinks=true)
    //-------------------------------------------------------------------------

    // ---- Run ---------------------------------------------------------------
    // generate index
    indexReferenceBWA(reference_fasta)
    // generate fasta file index
    generateFastaIndex(reference_fasta)
    // generate dict
    generateFastaDict(reference_fasta)

  emit:
    bwa_index_fls = indexReferenceBWA.out
    fasta_index_fl = generateFastaIndex.out
    dict_fl = generateFastaDict.out

}

/*
// -------------------------- DOCUMENTATION -----------------------------------
Indexes files are needed to do alignement with BWA MEM.
BWA requires .pac, .bwt, .ann, .amb and .sa index files that all have the same
basename (the full fasta name). Another file is an .alt index

Fasta Index are needed [TODO: add description]

Fasta Dict are needed because [TODO: add description]


https://gatk.broadinstitute.org/hc/en-us/articles/360037498992--How-to-Map-reads-to-a-reference-with-alternate-contigs-like-GRCH38
*/
