process genotype_amino_alleles {
    /*
    *                                           
    */

    input:
        path(vcf_file)

    output:   
        path("*.genotypes")

    script:
        genotypes_file = "${vcf_file.simpleName}.genotypes"
        """
        cat ${vcf_file} > ${genotypes_file}
        """
}
