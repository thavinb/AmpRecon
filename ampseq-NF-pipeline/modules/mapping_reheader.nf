process mapping_reheader {
     /*
     Run a Python script to copy across specific
     headers from the input bam file to the mapped bam file.
     */

    input:
        tuple val(sample_tag), path(selected_bam), path(selected_bam_metrics), path(clipped_bam)
        //val(sample_tag_clipped_bam) // sample tag provided by adapter clip
        //path(clipped_bam) // original cram with clipped adapters
        //val(sample_tag)
        //path(selected_bam) //new mapped bam
        path(reference_fasta) // reference fasta
        path(ref_dict)     // list of ref reference index files to stage

    output:
        val(sample_tag), emit: sample_tag
        path("${output_file}"), emit: reheadered_bam

    script:
        base_name=selected_bam.baseName
        output_file="${base_name}.reheadered.bam"
        """
        set -e
        set -o pipefail
        python3 ${projectDir}/modules/merge_headers.py ${clipped_bam} \
            ${selected_bam} \
            ${reference_fasta} | samtools reheader - \
            ${selected_bam} \
            | samtools merge ${output_file} -
        """
}
