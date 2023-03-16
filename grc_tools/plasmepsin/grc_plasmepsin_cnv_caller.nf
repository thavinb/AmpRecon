params.plasmepsin_cnv_output_filename = "plasmepsin_cnv_calls.tsv"

process grc_plasmepsin_cnv_caller {
    /*
    * Calls the breakpoint of the Plasmepsin 2/3 amplification in each of the supplied genotype files.
    * Uses plasmepsin genotypes, variants and loci supplied in the configuration file.
    * Writes these plasmepsin copy number variation calls to a single output tab-separated file.
    */

    input:
        val(genotype_file_list)

    output:
        path("${plasmepsin_out_filename}")

    script:
        grc_settings = params.grc_settings_file_path
        plasmepsin_out_filename = params.plasmepsin_cnv_output_filename
        genotype_file_list = genotype_file_list.join(" ")

        """
        grc_plasmepsin_cnv_caller.py \
            --genotype_files ${genotype_file_list} \
            --output_file "${plasmepsin_out_filename}" \
            --config "${grc_settings}"
        """
}
