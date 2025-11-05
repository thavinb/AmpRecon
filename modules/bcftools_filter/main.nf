// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.


process BCFTOOLS_FILTER {
    /*
    SNPs in the input BCF file are filtered and output as an uncompressed VCF file
    */

    tag "${meta.uuid}"
    label "bcftools"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.8--0' : 
        'biocontainers/bcftools:1.8--0' }"

    input:
        tuple val(meta), path(input_bcf)

    output:
        tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
        path "versions.yml"                             , emit: versions

    script:
        def prefix = "${meta.uuid}"
        intermediate_vcf="${prefix}.intermediate.vcf"

    // Had to escape backslash character in FORMAT line of script

        """
        # mark low depth
        bcftools filter \
        --mode + \
        --soft-filter LowDepth \
        --exclude FMT/DP\\<8 \
        --output-type v \
        < ${input_bcf} \
        > "${intermediate_vcf}"

        # mark low quality
        bcftools filter \
        --mode + \
        --soft-filter LowQual \
        --exclude 'QUAL<15 || MQ<20' \
        --output-type z \
        < "${intermediate_vcf}" \
        > "${prefix}.vcf.gz"

        # Indexing
        bcftools index -t "${prefix}.vcf.gz"

        # remove temporary file
        rm "${intermediate_vcf}"

        cat <<-EOF > versions.yml
        "${task.process}": 
            bcftools: "\$( bcftools --version | grep -E '^bcftools' | cut -f2 -d ' ')"
            htslib: "\$( bcftools --version | grep -E '^Using htslib' | cut -f 3 -d ' ')"
        EOF
        """
}
