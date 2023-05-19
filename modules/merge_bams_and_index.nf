
process merge_bams_and_index {
    /*
     * Merges multiple bam files.
     */

    input:
        tuple val(sample_id), file(bam_files)

    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
        
        """
        input_files=(${bam_files})
        if [ 1 == \${#input_files[@]} ]; then
            cp ${bam_files} ${sample_id}.bam
        else
            samtools merge ${sample_id}.bam ${bam_files}
	    samtools index ${sample_id}.bam
        fi
        """
}
