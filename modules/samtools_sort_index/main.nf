// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process SAMTOOLS_SORT_INDEX {
    /*
    * Sort and Index SAM file
    */
    tag "${meta.uuid}"
    label 'samtools'

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
