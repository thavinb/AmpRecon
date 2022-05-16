process variant_calling {
    /*
    *
    */

    input:
        path(bam_file)

    output:
        path("*.vcf")

    script:
        vcf_file = "${bam_file.simpleName}.vcf"
        """
        mv ./${bam_file} ./${vcf_file}
        """
}
