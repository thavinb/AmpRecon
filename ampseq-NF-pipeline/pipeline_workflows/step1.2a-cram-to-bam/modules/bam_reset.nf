params.bamreset_resetaux = 0
params.bamreset_level = 9
params.bamreset_verbose = 0


process bam_reset {
    /*
    * resets a BAM file to a pre-aligned state
    */
    errorStrategy 'ignore'

    input:
        val(sample_tag)
        path(collated_bam)
        // GAMBIARRA ---------------------------------------------------------
        // I had to find a way to be sure the reference fasta follow the
        // sample_cram. Easiest solution is just make it follow it all the way
        // to the process it is needed.
        //path(ref_fasta_file)
        //path(ref_bwa_index_fls)
        //path(ref_fasta_idx_fl)
        //path(ref_fasta_dct)
        // -------------------------------------------------------------------

    output:
        val("${sample_tag}"), emit: sample_tag
        path("${base_name}.reset.bam"), emit: reset_bam
        // GAMBIARRA ---------------------------------------------------------
        //path(ref_fasta_file), emit: ref_fasta_file
        //path(ref_bwa_index_fls), emit: ref_bwa_index_fls
        //path(ref_fasta_idx_fl), emit: ref_fasta_idx_fl
        //path(ref_fasta_dct), emit: ref_fasta_dct
        // -------------------------------------------------------------------


    script:
        base_name=collated_bam.baseName
        """
        bamreset \
            resetaux=${params.bamreset_resetaux} \
            level=${params.bamreset_level} \
            verbose=${params.bamreset_verbose} \
            < "${collated_bam}" \
            > "${base_name}.reset.bam"
        """
}
