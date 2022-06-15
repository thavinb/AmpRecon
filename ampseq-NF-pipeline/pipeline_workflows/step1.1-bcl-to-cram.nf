#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { basecalls_conversion } from '../modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from '../modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'

workflow bcl_to_cram {
    take:
        bcl_dir
        lane
        read_group
        study_name
        barcode
    main:
        // convert basecalls
        basecalls_conversion(bcl_dir, lane, read_group, study_name)

        // decode multiplexed bam file
        decode_multiplexed_bam(basecalls_conversion.out, barcode)

        // find adapter contamination in bam
        bam_find_adapter(decode_multiplexed_bam.out.output_file)

        // split bam by read group into cram
        bam_to_cram(bam_find_adapter.out)
        cram_ch = bam_to_cram.out

    emit:
        cram_ch
}
