// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process SCRAMBLE_CRAM_TO_BAM {
    /*
    Converts a cram to bam.
    */
    tag "${meta.uuid}"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staden_io_lib:1.15.1--hfc9290b_0' : 
        'biocontainers/staden_io_lib:1.15.1--hfc9290b_0' }"

    input:
        tuple val(meta), path(cram)
    output:
        tuple val(meta), path("*.bam"), emit: bam
        path "versions.yml"           , emit: versions

    script:
        def prefix = "${meta.uuid}"

        """
        scramble \
            -t 1 \
            -7 \
            -r "${meta.reference.fasta}" \
            -I cram \
            -O bam \
            < "${cram}" \
            > "${prefix}.bam"

        cat <<-EOF > versions.yml
        "${task.process}": 
            staden_io: \$(scramble -h | grep sCRAMble | sed 's/.*version //g')
        EOF
        """
}

