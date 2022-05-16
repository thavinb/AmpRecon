process variant_calling {
    /*
    *                                           
    */

    input:
        tuple path(bam_file)

    output:   
        tuple path("*.vcf")

    script:
        vcf_file = "${bam_file.simpleName}.vcf"
        """
        mv ${bam_file} ${vcf_file}
        """
}
