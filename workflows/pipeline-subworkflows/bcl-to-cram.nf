#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { basecalls_conversion } from '../../modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from '../../modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from '../../modules/bam_find_adapter.nf'
include { bam_to_cram } from '../../modules/bam_to_cram.nf'
include { rename_cram_fls } from '../../modules/rename_cram_fls.nf'

workflow BCL_TO_CRAM {
    take:
        pre_process_input_ch
        manifest
    main:
        // convert basecalls
        basecalls_conversion(pre_process_input_ch)
        decode_In_ch = pre_process_input_ch.join(basecalls_conversion.out)
        // decode multiplexed bam file
        decode_multiplexed_bam(decode_In_ch)
        find_adapter_In_ch = decode_multiplexed_bam.out

        // find adapter contamination in bam
        bam_find_adapter(find_adapter_In_ch)
        // split bam by read group into cram
        bam_to_cram(bam_find_adapter.out.run_id,
                    bam_find_adapter.out.bam_adapter_file,
                    bam_find_adapter.out.bam_metrics_file)

        //cram_ch = bam_to_cram.out

        // rename samples to samplesheet provided names
        // on this step the pattern for file names are set as
        // [run_id]_[lane]#[index]_[sample_name]
        // pannels names should be added at CRAM_TO_BAM 
        rename_cram_fls(bam_to_cram.out.run_id,
                        bam_to_cram.out.metrics_bam_file,
                        bam_to_cram.out.cram_fls,
                        params.lane,
                        manifest
                        )
        cram_ch = rename_cram_fls.out // tuple (run_id, cram_file)

        // add new sample tag to cram_ch
        cram_ch
          | flatten()
          | map { it -> tuple(it.simpleName, it, params.run_id) }
          | set {final_cram_ch} // tuple(sample_tag, cram_fl, run_id)
        //final_cram_ch.first().view()
    emit:
        final_cram_ch // tuple( sample_tag, cram_fl, run_id)
        //manifest_out
}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Mixed-sample MiSeq run sequencing data is first basecalled and then written to a BAM file.
[2] Resulting multiplexed BAM files are demultiplexed to separate reads based upon their sample of origin.
[3] BAM files are scanned for adapter sequence contamination, information about which is used to populate several auxiliary fields.
[4] BAM files are split by read group into CRAM files, which are emitted ready for further processing.

*/

