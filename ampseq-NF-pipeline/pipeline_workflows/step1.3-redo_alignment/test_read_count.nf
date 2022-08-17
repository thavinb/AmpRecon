#!/usr/bin/env nextflow
include { read_count_per_region_qc } from './modules/read_count_per_region.nf'

workflow {
        run_id="99999"
        cram_dir_ch="${params.cram_dir}"
        qc_run_ids_ch = Channel.from("GRC1", "GRC2", "Spec")
        qc_run_cnf_files_ch = Channel.from(file(params.grc1_qc_file), file(params.grc2_qc_file), file(params.spec_qc_file))
        read_count_per_region_qc(
        run_id,
        cram_dir_ch,
        qc_run_ids_ch,
        qc_run_cnf_files_ch
        )
}
