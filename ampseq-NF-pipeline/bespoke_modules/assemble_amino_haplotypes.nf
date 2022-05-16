process assemble_amino_haplotypes {
    /*
    *                                           
    */

    input:
        path(genotypes_file)

    output:   
        path("*.haplotypes")

    script:
        haplotypes_file = "${genotypes_file.simpleName}.haplotypes"
        """
        cat ${genotypes_file} > ${haplotypes_file}
        """
}
