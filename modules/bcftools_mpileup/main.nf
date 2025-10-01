// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.


process BCFTOOLS_MPILEUP {
    /*
    Creates an uncompressed BCF file containing calculated genotype likelihoods for every possible genomic position supported by the BAM
    */

    tag "${meta.uuid}"
    label "bcftools"
    label "process_low"

    input:
        tuple val(meta), path(bam), path(bai)
        /* tuple val(sample_tag), path(input_bam), path(input_bam_index), val(reference_file), val(reference_annotation_vcf) */

    output:
        tuple val(meta), path("*.bcf"), emit: bcf
        path "versions.yml"           , emit: versions

    script:
        def prefix = "${meta.uuid}"
        /* base_name = input_bam.baseName */
        /* output_bcf="${base_name}.bcf" */
        def min_bq = params.mpileup_min_bq ? params.mpileup_min_bq : 20
        def max_depth = params.mpileup_max_depth ? params.mpileup_max_depth : 50000

        """
        bcftools mpileup \
            --min-BQ ${min_bq} \
            --max-depth ${max_depth} \
            --annotate FORMAT/AD,FORMAT/DP \
            --targets-file "${meta.snps}" \
            --fasta-ref "${meta.reference.fasta}" \
            --output-type u \
            "${bam}" \
            > "${prefix}.bcf"
            
        cat <<-EOF > versions.yml
        "${task.process}": 
            bcftools: \$( bcftools --version | grep -E "^bcftools|^Using htslib" | tr '\\n' '\\t' | sed 's/bcftools //g'  )
        EOF
        """
}
