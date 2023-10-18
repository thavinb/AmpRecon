// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process bam_to_fastq {
    /*
    Converts BAM files to FASTQ.
    */

    input:
        tuple val(sample_tag), path(input_bam)

    output:
        tuple val(sample_tag), path("${base_name}.fastq")

    script:
        base_name=input_bam.baseName
        """
        bamtofastq \
            < "${input_bam}" \
            > "${base_name}.fastq"
        """
}
