params.samtools='samtools'
params.java='java'
params.picard='/bin/picard.jar'
params.gatk='GenomeAnalysisTK.jar'
params.gatk_print_reads_options=''

process bqsr {
    label 'genotyping'
    input:
        tuple val(sample_tag), path(bam_file), path(bam_index_file), val(reference_file), path(snp_list)

    output:
        tuple val(sample_tag), path("$recalibrated_bam_file")

    script:
        base_name=bam_file.getBaseName()
        gatk_recalibration_report="${base_name}.grp"
        recalibrated_bam_file="${base_name}.recalibrated.bam"

        java=params.java
        gatk=params.gatk
        picard=params.picard

        gatk_print_reads_gatk3_v2_memory=3500
        gatk_print_reads_options=params.gatk_print_reads_options
        jvm_args="-Xmx${gatk_print_reads_gatk3_v2_memory}m" // ???

        """
        # Base recalibration
        ${java} ${jvm_args} -jar ${gatk} -T BaseRecalibrator -R ${reference_file} -I ${bam_file} -o ${gatk_recalibration_report} --knownSites ${snp_list}

        # GATK print reads
        ${java} ${jvm_args} -jar ${gatk} -T PrintReads -R ${reference_file} -I ${bam_file} -o ${recalibrated_bam_file} -BQSR ${gatk_recalibration_report} ${gatk_print_reads_options}
        """
}
