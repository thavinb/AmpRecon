// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.


process BCFTOOLS_CALL {
    /*
    Calls SNPs from a BCF file containing all possible genotype likelihoods across genome
    */

    tag "${meta.uuid}"
    label "process_low"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.8--0' : 
        'biocontainers/bcftools:1.8--0' }"

    input:
        tuple val(meta), path(input_bcf)

    output:
        tuple val(meta), path("*.bcf"), emit: bcf
        path "versions.yml"           , emit: versions 

    script:
        def prefix = "${meta.uuid}"
        base_name = input_bcf.baseName
        output_bcf="${base_name}.bcftools_genotyped.bcf"
        def ploidy = "* * * * 2"

        """
        # Create Plody File
        echo "${ploidy}" > ploidy_file.ploidy

        bcftools call \
            --multiallelic-caller \
            --keep-alts \
            --skip-variants indels \
            --ploidy-file "ploidy_file.ploidy" \
            --output-type u \
            < "${input_bcf}" \
            > "${prefix}.genotyped.bcf"

        cat <<-EOF > versions.yml
        "${task.process}": 
            bcftools: \$( bcftools --version | grep -E "^bcftools|^Using htslib" | tr '\\n' '\\t' | sed 's/bcftools //g'  )
        EOF
        """
}

