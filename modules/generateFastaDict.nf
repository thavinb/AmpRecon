process generateFastaDict {
    /*
    * Indexes reference fasta file using bwa.
    */
    publishDir "${params.results_dir}/reference_files/", overwrite: true, mode: "copy"

    input:
        path(reference_fasta)

    output:
        path("${reference_fasta}.dict")

    script:
        """
        java -jar /app/picard.jar \
         CreateSequenceDictionary \
            REFERENCE=${reference_fasta} \
            OUTPUT=${reference_fasta}.dict
        """
}
