// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process SAMTOOLS_FLAGSTATS {
    /*
    * FLAGSTATS bam file
    */
    tag "${meta.uuid}"
    label 'samtools'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.flagstats"), emit: flagstats
        path "versions.yml"                 , emit: versions

    script:
        def prefix   = "${meta.uuid}"
        """
        samtools flagstats ${bam} > ${prefix}.flagstats

        cat <<-EOF > versions.yml
        "${task.process}":
            samtools: \$( samtools --version | grep -E "^samtools|^Using htslib" | tr '\\n' '\\t' | sed 's/samtools //g'  )
        EOF
        """
}
