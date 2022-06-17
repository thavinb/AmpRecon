params.bwa='bwa'
params.bwa_num_threads = 24
params.bwa_batch_input_bases = 100000000

process align_bam {
    /*
    * Map reads to reference
    */
    container ''

    input:
        val(tag)
        path(input_fastq)
        val(reference)        // reference fasta file
        path(ref_files)       // list of ref reference index files

    output:
        tuple val(tag), path("${output_file}.sam")

    script:
        bwa=params.bwa
        output_file=input_fastq.getBaseName()

        """
        bwa mem \
            -p \
            -Y \
            -K ${params.bwa_batch_input_bases} \
            -t ${params.bwa_num_threads} \
            "${reference}" \
            "${input_fastq}" \
            > "${output_file}.sam"
        """

}
