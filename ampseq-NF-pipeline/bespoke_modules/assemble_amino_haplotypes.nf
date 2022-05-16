process assemble_amino_haplotypes {
    /*
    *                                           
    */

    input:
        tuple path(genotypes_file)

    output:   
        tuple path("*.haplotypes")

    script:
        haplotypes_file = "${genotypes_file.simpleName}.haplotypes"
        """
        cat ${genotypes_file} > ${haplotypes_file}
        """
}
