params.grc_1_name = "GRC1.tsv"

process grc_assemble {
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    label "grc_tools"
    input:
        path(kelch13_mutation_calls)
        path(plasmepsin_cnv_calls)
        path(barcodes)
        path(species_calls)
        path(coi_estimates)
        path(drl_haplotypes)

    output:
        path("${output_grc_1}")

    script:
        output_grc_1 = params.grc_1_name

        """
        grc_assemble.py -grcs_in ${kelch13_mutation_calls} ${plasmepsin_cnv_calls} ${barcodes} ${species_calls} ${coi_estimates} ${drl_haplotypes} -grc_out_name ${output_grc_1}
        """
}
