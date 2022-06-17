process mapping_reheader {
     /*
     Run a Python script to copy across specific
     headers from the input bam file to the mapped bam file.
     */

    input:
        val(tag)
        path(alignment_input_file)
        path(mapped_input_file_to_reheader)
        val reference_file                  // fasta
        path ref_files                      // list of ref reference index files to stage

    output:
        tuple val(tag), path("${output_file}")

    script:
        base_name=mapped_input_file_to_reheader.baseName
        output_file="${base_name}.reheadered.bam"

        """
        set -e
        set -o pipefail
        merge_headers.py "${alignment_input_file}" "${mapped_input_file_to_reheader}" "${reference_file}" | samtools reheader - \
            "${mapped_input_file_to_reheader}" \
            > "${output_file}"
        """
}

