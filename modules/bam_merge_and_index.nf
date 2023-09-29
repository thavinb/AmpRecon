process bam_merge_and_index {
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
            cp ${bam_files} ${sample_key}.merged.bam
	    else
            samtools merge ${sample_key}.merged.bam ${bam_files}
	    samtools index ${sample_key}.merged.bam
        fi
        """
}
