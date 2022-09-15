params.bwa='bwa'
params.bwa_num_threads = 24
params.bwa_batch_input_bases = 100000000

process align_bam {
    /*
    * Map reads to reference
    */
    //publishDir "${params.results_dir}", overwrite: true
    label 'bwa'

    input:
        tuple val(sample_tag), path(fastq), path(reference_fasta), path(ref_bwa_index_fls), val(pannel_name)

    output:
        val("${sample_tag}"), emit: sample_tag
        path("${basename}_${pannel_name}.sam"), emit: sam_file

    script:
        bwa=params.bwa
        basename=fastq.simpleName
        """
        bwa mem \
            -p \
            -Y \
            -K ${params.bwa_batch_input_bases} \
            -t ${params.bwa_num_threads} \
            "${reference_fasta}" \
            "${fastq}" \
            > "${basename}_${pannel_name}.sam"
        """
}
// --- | WARNING | ------------------------------------------------------------
// currently no .alt is used, this means that reads that can map multiply are
// given low or zero MAPQ scores.
// the lack of this files is equivalent to add the -j option to the command
// disables the alt-handling.

// TODO: decide if alt file should be added or not
// ----------------------------------------------------------------------------

/*
  --- | DOCUMENTATION | -------------------------------------------------------
 https://manpages.ubuntu.com/manpages/bionic/man1/bwa.1.html
 we use -p to signal this is an interleaved paired-end fastq
 we use -Y to signal soft-clipping CIGAR operation for supplementary alignments

 Information regarding soft-clpping and what CIGAR string is, look here
 https://sites.google.com/site/bioinformaticsremarks/bioinfo/sam-bam-format/what-is-a-cigar
  ------------------------------------------------------------------------------
*/

