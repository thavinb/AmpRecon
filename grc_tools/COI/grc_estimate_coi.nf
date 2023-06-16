params.coi_output_filename = "coi.tsv"

process grc_estimate_coi {
    label "coi"
    input:
        path(barcodes_grc)

    output:
        path("${coi_out_filename}")

    script:
        grc_settings = params.grc_settings_file_path
        coi_out_filename = params.coi_output_filename
        input_suffix = "myRun" // this could be the run id (or any sequence batch identifier)
        
        """
        # write McCOIL in
        grc_process_mccoil_io.py -write_mccoil_in \
                --barcodes_files ${barcodes_grc} --config ${grc_settings} \
                --output_file ${input_suffix}.tsv

        # run McCOIL (if it is running inside a container, be sure the path set is '/app/THEREALMcCOIL/')
        Rscript ${params.mccoil_repopath}runMcCOIL.R -i ${input_suffix}.tsv --outPrefix ${input_suffix} \
            --totalRun ${params.mccoil_ntotal} --totalBurnIn ${params.mccoil_nburnin} --seed ${params.mccoil_seed} \
            --maxCOI ${params.mccoil_maxCOI} --M0 ${params.mccoil_m0} \
            --maxMissTol ${params.mccoil_maxMissTol} --e1 ${params.mccoil_e1} --e2 ${params.mccoil_e2}

        # write coi grc
        grc_process_mccoil_io.py -write_coi_grc \
            --mccoil_sum_file ${input_suffix}_summary.txt \
            --output_file ${coi_out_filename}
        
        # clean mccoil iterations log file
        rm ${input_suffix}
        """
}
