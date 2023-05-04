process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */
    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"
    input:
        tuple val(sample_tag), path(fastq), val(reference_fasta)

    output:
        tuple val(sample_tag), path(bam_file), path("${bam_file}.bai")

    script:
        bam_file="${sample_tag}.bam"

        """
        set -e
        set -o pipefail

        bwa mem -p -Y -K 100000000 -t 1 \
            "${reference_fasta}" \
            "${fastq}" | \
            scramble -0 -I sam -O bam | \
            samtools sort --threads 2 -o ${bam_file}.sorted.bam 
        
        samtools index ${bam_file}.sorted.bam

            mv ${bam_file}.sorted.bam ${bam_file}
            mv ${bam_file}.sorted.bam.bai ${bam_file}.bai
        """
}
