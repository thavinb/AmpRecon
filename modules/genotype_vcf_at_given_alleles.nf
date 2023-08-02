params.tabix='tabix'
params.gatk='/bin/GenomeAnalysisTK.jar'
params.bcftools='bcftools'
params.bgzip='bgzip'
params.suffix='genotyped_gatk.vcf'

process genotype_vcf_at_given_alleles {
    publishDir "${params.results_dir}/vcfs/", overwrite: true, mode: "copy"

    label 'vcf_parsing'

    input:
        tuple val(sample_tag), path(gvcf_fn), path(gvcf_fn_index), val(reference_file), path(snp_list)

    output:
        tuple val(sample_tag), path("${output_vcf_gz}"), path(output_vcf_gz_index), emit: vcf_file_and_index

    script:        
        java_memory = task.memory.mega - 200

        tabix=params.tabix
        gatk=params.gatk
        bcftools=params.bcftools
        bgzip=params.bgzip
        suffix=params.suffix

        output_base = gvcf_fn.getBaseName().replace(".recalibrated","").replace(".vcf", "")
        genotyped_gvcf = "${output_base}.GenotypeGVCFs.vcf"

        output_vcf="${output_base}.${suffix}"
        output_vcf_gz="${output_vcf}.gz"
        output_vcf_gz_index="${output_vcf_gz}.tbi"

        """
        java -Xmx${java_memory}m \
        -jar ${gatk} \
        -T GenotypeGVCFs \
        -R ${reference_file} \
        --variant ${gvcf_fn} \
        -allSites \
        -L ${snp_list} \
        -o ${genotyped_gvcf}

        main ${genotyped_gvcf} \
        ${snp_list} \
        --output_suffix ${suffix}

        # Compress and index (g)vcf
        ${bgzip} -c ${output_vcf} > ${output_vcf_gz}
        ${tabix} -p vcf ${output_vcf_gz}

        rm ${genotyped_gvcf}*
        """
}
