params.samtools='samtools'
params.java='java'
params.picard='/bin/picard.jar'
params.gatk='GenomeAnalysisTK.jar'
params.gatk_print_reads_options=''
params.gatk_base_recalibrator_options=''

process bqsr {
    //label 'pf7_container'
    label 'genotyping'
    input:
        tuple val(sample_tag), path(bam_file), path(bam_index_file), path(reference_file), path(reference_index_file), path(reference_dict_file)

    output:
        tuple val(sample_tag), path("$recalibrated_bam_file")

    script:
        base_name_ref=reference_file.getBaseName()
        dict_file="${base_name_ref}.dict"
        base_name=bam_file.getBaseName()
        gatk_recalibration_report="${base_name}.grp"
        recalibrated_bam_file="${base_name}.recalibrated.bam"

        java=params.java
        gatk=params.gatk
        picard=params.picard

        gatk_print_reads_gatk3_v2_memory=3500
        gatk_print_reads_options=params.gatk_print_reads_options
        gatk_base_recalibrator_options=params.gatk_base_recalibrator_options
        jvm_args="-Xmx${gatk_print_reads_gatk3_v2_memory}m" // ???

        """
        mv ${reference_dict_file} ${dict_file}

        # Base recalibration
        ${java} ${jvm_args} -jar ${gatk} -T BaseRecalibrator -R ${reference_file} -I ${bam_file} -o ${gatk_recalibration_report} ${gatk_base_recalibrator_options}

        # GATK print reads
        ${java} ${jvm_args} -jar ${gatk} -T PrintReads -R ${reference_file} -I ${bam_file} -o ${recalibrated_bam_file} -BQSR ${gatk_recalibration_report} ${gatk_print_reads_options}
        """
}
