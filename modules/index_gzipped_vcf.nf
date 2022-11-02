process index_gzipped_vcf {
    publishDir "${params.results_dir}/", overwrite: true
    input:
        tuple val(sample_tag), path(input_zipped_vcf_file)

    output:
        tuple val(sample_tag), path(input_zipped_vcf_file), path("${vcf_name}.tbi")

    script:
        base_name = input_zipped_vcf_file.baseName
        vcf_name="${base_name}.gz"

        """
           bcftools index --tbi "${vcf_name}" 
        """
}
