params.bwa='bwa'
params.bwa_num_threads = 1
params.bwa_batch_input_bases = 100000000

process bwa_alignment {
    /*
    * Map reads to reference
    */
    //publishDir "${params.results_dir}", overwrite: true
    label 'bwa'

    input:
        tuple val(sample_tag), path(fastq), val(reference_fasta), val(panel_name)

    output:
        val("${sample_tag}"), emit: sample_tag
        path("${sample_tag}-${ref_simplename}.sam"), emit: sam_file

    script:
        bwa=params.bwa
        ref_simplename=file(reference_fasta).simpleName
        """
        bwa mem \
            -p \
            -Y \
            -K ${params.bwa_batch_input_bases} \
            -t ${params.bwa_num_threads} \
            "${reference_fasta}" \
            "${fastq}" \
            > "${sample_tag}-${ref_simplename}.sam"
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

