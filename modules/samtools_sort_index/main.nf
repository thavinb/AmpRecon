// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process SAMTOOLS_SORT_INDEX {
    /*
    * Sort and Index SAM file
    */
    tag "${meta.uuid}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.6--h5fe306e_13' : 
        'biocontainers/samtools:1.6--h5fe306e_13' }"

    input:
        tuple val(meta), path(sam)

    output:
        tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
        path "versions.yml"                          , emit: versions

    script:
        def prefix   = "${meta.uuid}"
        """
        set -e
        set -o pipefail

        samtools view \\
            -bu \\
            --threads ${task.cpus} \\
            "${sam}" | \\
            samtools sort -n - | \\
            samtools fixmate - - | \\
            samtools sort -o ${prefix}.bam 

        samtools index ${prefix}.bam

        cat <<-EOF > versions.yml
        "${task.process}":
            samtools: \$( samtools --version | grep -E "^samtools|^Using htslib" | tr '\\n' '\\t' | sed 's/samtools //g'  )
        EOF
        """
}
