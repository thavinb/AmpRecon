params.samtools = 'samtools'
params.bgzip = 'bgzip'
params.gatk = 'gatk'
params.tabix = 'tabix'
params.java='java'

params.gatk_haplotype_caller_gatk4_contamination = 0
params.gatk_haplotype_caller_gatk4_interval_file = ''

process gatk_haplotype_caller_gatk4 {

    label 'genotyping'
    input:
        tuple val(sample_tag), path(bam_file), val(reference_file)

    output:
        tuple val(sample_tag), path("${vcf_file_gz}"), path("${vcf_file_gz_index}"), emit: vcf_file_and_index

    script:
        java_Xmx = 16000
        jvm_args="-Xmx${java_Xmx}m"

        base_name=bam_file.getBaseName()
        vcf_file="${base_name}.vcf"
        vcf_file_gz="${vcf_file}.gz"
        vcf_file_gz_index="${vcf_file_gz}.tbi"

        samtools = params.samtools
        bgzip = params.bgzip
        java = params.java
        gatk = params.gatk
        tabix = params.tabix
        contamination = params.gatk_haplotype_caller_gatk4_contamination

        if (params.gatk_haplotype_caller_gatk4_interval_file != '') {
            interval_file = '-L ${params.gatk_haplotype_caller_gatk4_interval_file}'
        }
        else {
            interval_file = ''
        }

        """
        ${samtools} index ${bam_file}

        # GATK4 Haplotyper has issues with OpenJDK 11

        # Run haplotype caller
        ${gatk} \
        --java-options ${jvm_args} \
        HaplotypeCaller \
        -R ${reference_file} \
        -I ${bam_file} \
        -O ${vcf_file} \
        -contamination ${contamination} \
        ${interval_file} \
        -ERC GVCF

        # Compress and index (g)vcf
        ${bgzip} -c ${vcf_file} > ${vcf_file_gz}
        ${tabix} -p vcf ${vcf_file_gz}
        """
}
