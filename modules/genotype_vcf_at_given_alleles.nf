params.tabix='tabix'
//params.alleles_fn='/lustre/scratch118/malaria/team112/pipelines/resources/pfalciparum/pf_6x_genotyping.vcf.gz'
params.vcf_format='{sample}.GenotypeGVCFs.{alleles}.vcf.gz'
params.gatk='GenomeAnalysisTK.jar'
params.bcftools='bcftools'
params.bgzip='bgzip'

process genotype_vcf_at_given_alleles {
    /**
    * Runs the genotype_vcf_at_given_alleles.py script from docker
    * Example call from https://notebooks.gesis.org/binder/jupyter/user/ipython-ipython-in-depth-47gm2yay/notebooks/binder/20191112_genotype_pf$

    ./genotype_vcf_at_given_alleles.py \
    --gvcf_fn /lustre/scratch116/malaria/pipelines/setups/pf_63/output/1/6/a/4/1363815/1_gatk_haplotype_caller_gatk3_v3/gatk.vcf.gz \
    --alleles_fn /lustre/scratch118/malaria/team112/pipelines/resources/pfalciparum/pf_6x_genotyping.vcf.gz \
    --genotyped_fn test_RCN11295.pf_6x.vcf \
    --GenotypeGVCFs_format test_{sample}.GenotypeGVCFs.{alleles}.vcf.gz \
    --overwrite

    */
    //label 'pf7_container'
    label 'genotyping'
    
    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/genotype_vcf_at_given_alleles:0.0.3'

    input:
        tuple val(sample_tag), path(gvcf_fn), path(gvcf_fn_index), path(reference_file), path(reference_index_file), path(reference_dict_file), path(snp_list)

    output:
        tuple val(sample_tag), path("${output_vcf_gz}"), emit: vcf_file

    script:
        base_name_ref=gvcf_fn.getBaseName().replace(".recalibrated.vcf", "")
        output_vcf="${base_name_ref}.gatk_genotyped.vcf"
        output_vcf_gz="${output_vcf}.gz"
        output_vcf_gz_index="${output_vcf_gz}.tbi"
        java_memory = 8000 * task.attempt

//        alleles_fn=params.alleles_fn
        vcf_format=params.vcf_format
        tabix=params.tabix
        gatk=params.gatk
        bcftools=params.bcftools
        bgzip=params.bgzip

        base_name_ref=reference_file.getBaseName()
        dict_file="${base_name_ref}.dict"

        """
        if [["${reference_dict_file}" != "${dict_file}"]]; then mv ${reference_dict_file} ${dict_file}; fi

        genotype_vcf_at_given_alleles.py \
        --java_memory ${java_memory}m \
        --gvcf_fn ${gvcf_fn} \
        --alleles_fn ${snp_list} \
        --genotyped_fn ${output_vcf} \
        --GenotypeGVCFs_format ${vcf_format} \
        --ref_fasta ${reference_file} \
        --gatk_jar ${gatk} \
        --bcftools ${bcftools} \
        --bgzip ${bgzip} \
        --overwrite
        """
}
