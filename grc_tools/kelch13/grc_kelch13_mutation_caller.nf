params.kelch13_mutations_output_filename = "kelch13_mutation_calls.tsv"

process grc_kelch13_mutation_caller {
    /*
    * Calls non-synonymous mutations in asupplied region of the kelch13 gene in each of the supplied genotype files.
    * Uses a file detailing the kelch13 reference sequence, codon structure and amino acid translation.
    * It uses a codon key file for translating nucleotide codons into associated amino acids
    * It also uses the kelch13 region supplied in the configuration file.
    * Writes these kelch13 mutation calls to a single output tab-separated file.
    */

    input:
        val(genotype_file_list)
        path(kelch_reference_file)
        path(codon_key_file)

    output:
        path("${kelch13_out_filename}")

    script:
        grc_settings = params.grc_settings_file_path
        kelch13_out_filename = params.kelch13_mutations_output_filename
        genotype_file_list = genotype_file_list.join(" ")

        """
        grc_kelch13_mutation_caller.py \
            --genotype_files ${genotype_file_list} \
            --output_file "${kelch13_out_filename}" \
            --config "${grc_settings}" \
            --kelch_reference_file "${kelch_reference_file}" \
            --codon_key_file "${codon_key_file}"
        """
}
