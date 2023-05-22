process merge_bams_and_index {
    /*
     * Merges multiple bam files.
     */

    input:
        tuple val(sample_key), file(bam_files)

    output:
        tuple val(sample_key), path("${sample_key}.bam"), path("${sample_key}.bam.bai")

    script:
        
        """
        input_files=(${bam_files})
        if [ 1 == \${#input_files[@]} ]; then
            cp ${bam_files} ${sample_key}.bam
	    samtools index ${sample_key}.bam 
        else
            samtools merge ${sample_key}.bam ${bam_files}
	    samtools index ${sample_key}.bam
        fi
        """
}
