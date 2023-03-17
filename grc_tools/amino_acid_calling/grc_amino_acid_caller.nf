params.drl_haplotypes_output_filename = "drl_haplotype_calls.tsv"
params.grc_2_filename = "GRC2.tsv"

process grc_amino_acid_caller {
    publishDir "${params.results_dir}/", pattern: "${output_grc_2}", overwrite: true, mode: "copy"
    label "grc_tools"
    input:
        val(genotype_files)
        path(drl_information_file)
        path(codon_key_file)

    output:
        path("${output_drl_haplotypes}"), emit: drl_haplotypes
        path("${output_grc_2}"), emit: grc2

    script:
        grc_settings = params.grc_settings_file_path
        output_drl_haplotypes = params.drl_haplotypes_output_filename
        output_grc_2 = params.grc_2_filename
        genotype_files = genotype_files.join(" ")

        """
        grc_amino_acid_caller.py \
            --genotype_files ${genotype_files} \
            --config "${grc_settings}"  \
            --output_grc1_file "${output_drl_haplotypes}" \
            --output_grc2_file "${output_grc_2}" \
            --drl_information_file "${drl_information_file}"  \
            --codon_key_file "${codon_key_file}"
        """

}
