process mapping_reheader {
     /*
     Run a Python script to copy across specific
     headers from the input bam file to the mapped bam file.
     */
    label 'python_plus_samtools'
    input:
        tuple val(sample_tag), path(scrambled_bam),  path(clipped_bam), path(reference_fasta), path(ref_dict)

    output:
        tuple val(sample_tag), path("${output_file}")

    script:
        base_name_ref=reference_fasta.getBaseName()
        dict_file="${base_name_ref}.fasta.dict"

        base_name=scrambled_bam.simpleName
        output_file="${sample_tag}.reheadered.bam"
        """
        # If needed, rename dictionary file to correct format.
        if [[ "${ref_dict}" != "${dict_file}" ]]; then mv "${ref_dict}" "${dict_file}"; fi

        set -e
        set -o pipefail
        merge_headers.py ${clipped_bam} \
            ${scrambled_bam} \
            ${reference_fasta} | samtools reheader - \
            ${scrambled_bam} \
            | samtools merge ${output_file} -
        """
}

