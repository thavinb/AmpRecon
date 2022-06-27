process generateFastaIndex {
    /*
    * Indexes reference fasta file using bwa.
    */
    publishDir "${params.results_dir}/reference_files/", overwrite: true

    input:
        path(reference_fasta)

    output:
        path("${reference_fasta}*")

    script:
        """
        samtools faidx ${reference_fasta}
        """
}
