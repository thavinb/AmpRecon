#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include { basecalls_conversion } from '../modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from '../modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'
include { rename_cram_fls } from '../modules/rename_cram_fls.nf'

include { collate_alignments } from '../modules/collate_alignments.nf'
include { bam_reset } from '../modules/bam_reset.nf'
include { clip_adapters } from '../modules/clip_adapters.nf'
include { bam_to_fastq } from '../modules/bam_to_fastq.nf'
include { bwa_alignment } from '../modules/bwa_alignment.nf'
include { scramble_sam_to_bam } from '../modules/scramble.nf'
include { mapping_reheader } from '../modules/mapping_reheader.nf'
include { bam_split } from '../modules/bam_split.nf'
include { bam_merge } from '../modules/bam_merge.nf'
include { alignment_filter } from '../modules/alignment_filter.nf'
include { sort_bam } from '../modules/sort_bam.nf'

// - process to extract and validate information expected based on input params
include { get_taglist_file } from '../modules/manifest2tag.nf'
include { make_samplesheet_manifest } from '../modules/make_samplesheet_manifest.nf'
include { validate_samplesheet_manifest } from '../modules/samplesheet_manifest_validation.nf'
include { PARSE_PANEL_SETTINGS } from './parse_panels_settings.nf'
include { retrieve_miseq_run_from_s3 } from '../modules/retrieve_miseq_run_from_s3.nf'

workflow BCL_TO_COLLATED_CRAM {
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

        // rename samples to samplesheet provided names
        // on this step the pattern for file names are set as
        // [run_id]_[lane]#[index]_[sample_name]
        // panels names should be added at CRAM_TO_BAM 
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
          | set {final_cram_ch} // tuple(file_id, cram_fl, run_id)

    emit:
        final_cram_ch // tuple( file_id, cram_fl, run_id)

}

workflow COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ {
    take:
        // manifest from step 1.1
        //intermediate_csv
        cram_ch // tuple(file_id, cram_fl, run_id) | file_id = [run_id]_[lane]#[index]_[sample_name]-
        file_id_reference_files_ch // tuple (sample_id, ref_fasta, panel_name)

    main:
        // --| DEBUG |--------------------------------------------------
        if (!(params.DEBUG_takes_n_bams == null)){
            collate_alignments(cram_ch.take(params.DEBUG_takes_n_bams))
        }
        // -------------------------------------------------------------
        else {
            collate_alignments(cram_ch)
        }

        // Transform BAM file to pre-aligned state
        bam_reset(collate_alignments.out.sample_tag,
                  collate_alignments.out.collated_bam)

        // Remove adapters
        clip_adapters(bam_reset.out.sample_tag,
                 bam_reset.out.reset_bam)

        // Convert BAM to FASTQ
        bam_to_fasq_In_ch = clip_adapters.out
                              | multiMap {it -> file_id: it[0]
                                                bam_clipped_file: it[1]
                                        }
        bam_to_fastq(bam_to_fasq_In_ch.file_id,
                    bam_to_fasq_In_ch.bam_clipped_file)

        // Align reads to reference genome
        bam_to_fastq.out //tuple (file_id, fastq_files)
              // get panel resource files
              | join(file_id_reference_files_ch) // tuple (file_id, fastq_files, panel_name, ref_fasta)
              | map{it -> tuple(it[0], it[1], it[2], it[3])} // tuple (new_file_id, fastq, fasta, panel_name)
              | set{ bwa_ch }

        bwa_alignment(bwa_ch)

        // Convert SAM to BAM
     
        scramble_sam_to_bam(bwa_alignment.out.sample_tag,
                            bwa_alignment.out.sam_file,
        )
        
        // Merges the current headers with the old ones.
        // Keeping @SQ.*\tSN:([^\t]+) matching lines from the new header.
        reheader_in_ch = scramble_sam_to_bam.out 
                            | join(clip_adapters.out)
                            | join(file_id_reference_files_ch)
                            | map { it -> tuple(it[0], it[1], it[2], it[3]) } // tuple(file_id, scrambled_bam, clipped_bam, ref_fasta)
        mapping_reheader(reheader_in_ch)

        // Split BAM rank pairs to single ranks per read
        bam_split(mapping_reheader.out)

        // Merge BAM files with same reads

        bam_merge_In_ch = bam_split.out.join(clip_adapters.out)

        bam_merge(bam_merge_In_ch)

        // Split alignments into different files
        alignment_filter(bam_merge.out.sample_tag,
                         bam_merge.out.merged_bam)

        // BAM sort by coordinate
        sort_bam(params.run_id,
                 alignment_filter.out.sample_tag,
                 alignment_filter.out.selected_bam)
        bam_ch = sort_bam.out
  
    emit:
        bam_ch
}

workflow MISEQ_TO_READS {
  take:
    reference_ch // tuple (fasta, panel_name, snp_list)

  main:
      
    if (params.s3_bucket_input != null) {
      retrieve_miseq_run_from_s3(params.s3_uuid)
      input_csv_ch = retrieve_miseq_run_from_s3.out.tuple_ch
    }
  
    else {
	   //create input channel if not running s3 entrypoint
      input_csv_ch = Channel.of(tuple(params.run_id,
        params.bcl_dir,
        params.lane,
        params.study_name,
        params.read_group,
        params.library))
    }
    // process samplesheets manifest (necessary to get barcodes) and validate it
    input_csv_ch
      | map {it -> tuple (it[0], it[1])} // tuple(run_id, bcl_dir)
      | set { make_samplesheet_In_ch}

    make_samplesheet_manifest(make_samplesheet_In_ch)
    panel_names_list = reference_ch.map{it -> it[1].toString()}.collect()

    validate_samplesheet_manifest(make_samplesheet_manifest.out.tuple, panel_names_list)

    get_taglist_file_In_ch = input_csv_ch.join(validate_samplesheet_manifest.out)
    get_taglist_file(get_taglist_file_In_ch)

    step1_Input_ch = input_csv_ch.join(get_taglist_file.out)

    // Stage 1 - Step 1: BCL to CRAM
    BCL_TO_COLLATED_CRAM(step1_Input_ch, make_samplesheet_manifest.out.manifest_file)
    cram_ch = BCL_TO_COLLATED_CRAM.out // tuple (file_id, cram_fl, run_id)

    // get the relevant sample data from the manifest
    file_id_ch = make_samplesheet_manifest.out.manifest_file // WARN: this need to be removed, we should no rely on results dir
      | splitCsv(header: ["lims_id", "sims_id", "index", "assay",
          "barcode_sequence", "well", "plate"],
          skip: 18)
      | map { row -> tuple(row.lims_id, row.assay, row.index) }
      | map{it -> tuple("${params.run_id}_${params.lane}#${it[2]}_${it[0]}", it[1], it[0])} // tuple (file_id, panel_name, lims_id)

    // link file_id to sample_id
    file_id_ch.map{it -> tuple("${it[0]}_${it[1]}", it[2])}.set{file_id_to_sample_id_ch}

    // add panel names to file_ids
    panel_name_cram_ch = cram_ch.join(file_id_ch) // tuple (file_id, cram_fl, run_id, panel_name)
     | map{it -> tuple("${it[0]}_${it[3]}", it[1], it[2])} // tuple(file_id_panel_name, cram_fl, run_id)
    new_file_id_ch = file_id_ch.map{it -> tuple("${it[0]}_${it[1]}", it[1])} // tuple (file_id, panel_name)  

    // assign each sample tag the appropriate set of reference files
    new_file_id_ch
      | combine(reference_ch,  by: 1) // tuple (panel_name, file_id, fasta,snp_list)
      | map{it -> tuple(it[1], it[0], it[2], it[3])}
      | set{file_id_reference_files_ch}
    // tuple(file_id, panel_name path/to/reference/genome, snp_list)

    // Stage 1 - Step 2: CRAM to BAM
    COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ(panel_name_cram_ch, file_id_reference_files_ch.map{it -> tuple(it[0], it[2], it[1])})  // tuple (new_file_id, ref_fasta, panel_name)
    bam_files_ch = COLLATED_CRAM_TO_SPLIT_CRAM_AND_FASTQ.out.bam_ch

  emit:
    bam_files_ch // tuple (file_id, bam_file)
    file_id_reference_files_ch // tuple (file_id, panel_name, path/to/reference/genome, snp_list)
    file_id_to_sample_id_ch // tuple (file_id, sample_id)
}


def count_miseq_to_reads_params_errors(){
    /*
    This functions counts the number of errors on input parameters exclusively used on Incountry workflow
    
    checks:
     - if irods manifest was provided
     - if irods manifest provided exists
    
    False for any of those conditions counts as an error.
    
    Returns
    ---
    <int> the number of errors found
    */

    def err = 0

    if (params.run_id == null){
      log.error("A run_id parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }

    if (params.bcl_dir == null && params.s3_bucket_input == null){
      log.error("Either a bcl directory or a s3 bucket input must be specified for in-country execution_mode.")
      err += 1
    } 
    

    if (params.lane == null){
      log.error("A lane parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    
    if (params.study_name == null){
      log.error("A study_name parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    if (params.read_group == null){
      log.error("A read_group parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    
    if (params.library == null){
      log.error("A library parameter must be provided for execution mode '${params.execution_mode}'.")
      err += 1
    }
    return err
}
