process mapping_reheader {
     /*
     Run a Python script to copy across specific
     headers from the input bam file to the mapped bam file.
     */

    input:
        tuple val(sample_tag), path(scrambled_bam),  path(clipped_bam)
        //path(clipped_bam) // original cram with clipped adapters
        //path(scrambled_bam) //new mapped bam
        path(reference_fasta) // reference fasta
        path(ref_dict)     // list of ref reference index files to stage

    output:
        tuple val(sample_tag), path("${output_file}")

    script:
        base_name=scrambled_bam.baseName
        output_file="${base_name}.reheadered.bam"
        """
        set -e
        set -o pipefail
        python3 ${projectDir}/modules/merge_headers.py ${clipped_bam} \
            ${scrambled_bam} \
            ${reference_fasta} | samtools reheader - \
            ${scrambled_bam} \
            | samtools merge ${output_file} -
        """
}
