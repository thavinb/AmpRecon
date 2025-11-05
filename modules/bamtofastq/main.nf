// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process BAMTOFASTQ {
    /*
    Converts BAM files to FASTQ.
    */

    tag "${meta.uuid}"
    label 'biobambam'
    label 'process_low'

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.fastq"), emit: fastq
        path "versions.yml"             , emit: versions

    script:
        def prefix = "${meta.uuid}"
        """
        bamtofastq \
            filename="${bam}" \
            F="${prefix}_1.fastq" \
            F2="${prefix}_2.fastq"

        cat <<-EOF > versions.yml
        "${task.process}":
            biobambam: \$(bamtofastq -v 2>&1 | head -n1 | sed -n 's/.*version //p' | sed 's/\\.\$//')
        EOF
        """
}
