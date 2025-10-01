// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.


process BCFTOOLS_STATS {
    /*
    Calls SNPs from a BCF file containing all possible genotype likelihoods across genome
    */

    tag "${meta.uuid}"
    label "bcftools"
    label "process_low"

    input:
        tuple val(meta), path(vcf), path(tbi)

    output:
        tuple val(meta), path("*.stats"), emit: stats
        path "versions.yml"             , emit: versions 

    script:
        def prefix = "${meta.uuid}"

        """
        bcftools stats ${vcf} > ${prefix}.stats

        cat <<-EOF > versions.yml
        "${task.process}": 
            bcftools: \$( bcftools --version | grep -E "^bcftools|^Using htslib" | tr '\\n' '\\t' | sed 's/bcftools //g'  )
        EOF
        """
}

