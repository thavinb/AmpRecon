process mapping_reheader {
     /*
     Run a Python script to copy across specific
     headers from the input bam file to the mapped bam file.
     */
    label 'python_plus_samtools'
    input:
        tuple val(sample_tag), path(scrambled_bam),  path(clipped_bam), val(reference_fasta)

    output:
        tuple val(sample_tag), path("${output_file}")

    script:

        base_name=scrambled_bam.simpleName
        output_file="${sample_tag}.reheadered.bam"
        """
        set -e
        set -o pipefail
        merge_headers.py ${clipped_bam} \
            ${scrambled_bam} \
            ${reference_fasta} | samtools reheader - \
            ${scrambled_bam} \
            | samtools merge ${output_file} -
        """
}

