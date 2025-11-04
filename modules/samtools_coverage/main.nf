// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process SAMTOOLS_COVERAGE {
    /*
    * Get coverage from bam file
    */
    tag "${meta.uuid}"
    label 'samtools'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.coverage"), emit: coverage
        path "versions.yml"                , emit: versions

    script:
        def prefix   = "${meta.uuid}"
        """
        samtools coverage ${bam} > ${prefix}.coverage

        cat <<-EOF > versions.yml
        "${task.process}":
            samtools: "\$( samtools --version | grep -E "^samtools" | cut -f2 -d ' '  )"
            htslib: "\$( samtools --version | grep -E '^Using htslib' | cut -f 3 -d ' ')"
        EOF
        """
}
