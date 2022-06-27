params.bwa='bwa'

process indexReferenceBWA {
    /**
    * Indexes reference fasta file using bwa.
    */
    publishDir "${params.results_dir}/reference_files/", overwrite: true

    input:
        path(reference_fasta)

    output:
        path("${reference_fasta}*")

    script:
        bwa=params.bwa

        """
        ${bwa} index -a bwtsw -p ${reference_fasta} ${reference_fasta}
        """
}
