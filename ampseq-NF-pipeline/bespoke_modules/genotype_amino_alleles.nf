process genotype_amino_alleles {
    /*
    *                                           
    */

    input:
        tuple path(vcf_file)

    output:   
        tuple path("*.genotypes")

    script:
        genotypes_file = "${vcf_file.simpleName}.genotypes"
        """
        cat ${vcf_file} > ${genotypes_file}
        """
}
