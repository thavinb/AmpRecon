#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
//include { ALIGN_BAM as ALIGNMENT } from '../../modules/align_bam.nf'
include { BAMBI_I2B } from '../../modules/bambi_i2b.nf'
include { BAM_TO_FASTQ } from '../../modules/bam_to_fastq.nf'
include { BAMCOLLATE2 } from '../../modules/bamcollate2.nf'
include { BAMRESET } from '../../modules/bam_reset.nf'
include { CLIP_ADAPTERS } from '../../modules/clip_adapters.nf'
include { GENERATE_FASTA_INDEX } from '../../modules/generateFastaIndex.nf'
include { SAMTOOLS_SPLIT } from '../../modules/samtools_split.nf'
include { BAMBI_DECODE } from '../../modules/bambi_decode.nf'
include { BAMADAPTERFIND } from '../../modules/bamadapterfind.nf'
include { RENAME_CRAM_FILES } from '../../modules/rename_cram_fls.nf'
include { SCRAMBLE_SAM_TO_BAM } from '../../modules/scramble.nf'

workflow CORE_PIPELINE_REPLICA {
    take:
        input_ch // extracted from input csv from user in EXTRACT_PARAMS
        barcodes_file
        tag_list_file
        reference_information

    main:
        // convert basecalls
        BAMBI_I2B(input_ch)

        BAMBI_DECODE(BAMBI_I2B.out, barcodes_file)
        bam_adapter_ch = BAMBI_DECODE.out.map { it -> tuple ( it[0], it[1] ) }

        // find adapter contamination in bam
        BAMADAPTERFIND(bam_adapter_ch)
        // split bam by read group into cram
        SAMTOOLS_SPLIT(BAMADAPTERFIND.out)

        RENAME_CRAM_FILES(SAMTOOLS_SPLIT.out)
        cram_ch = RENAME_CRAM_FILES.out

        BAMCOLLATE2(cram_ch)
        BAMRESET(BAMCOLLATE2.out)
        CLIP_ADAPTERS(BAMRESET.out)

        fastq_ch = BAM_TO_FASTQ(CLIP_ADAPTERS.out)

        //generate reference index
        reference_fasta = reference_information.map { it -> it[1]}
        reference_idx = GENERATE_FASTA_INDEX(reference_information)

        pre_alignment_ch = fastq_ch.join(reference_information)
        pre_alignment_ch = pre_alignment_ch.join(reference_idx)

    emit:
      cram_ch



}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Mixed-sample MiSeq run sequencing data is first basecalled and then written to a BAM file.
[2] Resulting multiplexed BAM files are demultiplexed to separate reads based upon their sample of origin.
[3] BAM files are scanned for adapter sequence contamination, information about which is used to populate several auxiliary fields.
[4] BAM files are split by read group into CRAM files, which are emitted ready for further processing.

*/
