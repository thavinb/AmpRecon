params.bwa='bwa'
params.bwa_num_threads = 24
params.bwa_batch_input_bases = 100000000

process align_bam {
    /*
    * Map reads to reference
    */
    //publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        val(sample_tag)
        path(fastq)
        path(reference_fasta) // fasta file
        path(ref_bwa_index_fls) // index files for the reference
        //val(tag)
        //path(input_fastq)
        //val(reference)        // reference fasta file
        //path(ref_files)       // list of ref reference index files

    output:
        val("${sample_tag}"), emit: sample_tag
        path("${basename}.sam"), emit: sam_file

    script:
        bwa=params.bwa
        basename=fastq.baseName
        """
        bwa mem \
            -p \
            -Y \
            -K ${params.bwa_batch_input_bases} \
            -t ${params.bwa_num_threads} \
            "${reference_fasta}" \
            "${fastq}" \
            > "${basename}.sam"
        """
}
// --- | WARNING | ------------------------------------------------------------
// currently no .alt is used, this means that reads that can map multiply are
// given low or zero MAPQ scores.
// the lack of this files is equivalent to add the -j option to the command
// disables the alt-handling.

// TODO: decide if alt file should be added or not
// ----------------------------------------------------------------------------
