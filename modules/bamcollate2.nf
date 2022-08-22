params.bamcollate_collate = 1
params.bamcollate_level = 9

process BAMCOLLATE2 {
    /*
    * collates bam/cram reads/alignments by name
    */
    errorStrategy 'ignore'

    //publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        //tuple val(run_id), path(cram_fl), val(sample_tag)
        tuple val(run_id), path(sample_cram)
        // GAMBIARRA ---------------------------------------------------------
        // I had to find a way to be sure the reference fasta follow the
        // sample_cram. Easiest solution is just make it follow it all the way
        // to the process it is needed.
        //path(ref_fasta_file)
        //path(ref_bwa_index_fls)
        //path(ref_fasta_idx_fl)
        //path(ref_fasta_fasta_dict)
        // -------------------------------------------------------------------
        //val(tag) // be sure what tag supose to mean here (not run id a supose)
        //path(input_file)

    output:
        tuple val(run_id), path("${base_name}.collated.bam")
        // GAMBIARRA ---------------------------------------------------------
        //path(ref_fasta_file), emit: ref_fasta_file
        //path(ref_bwa_index_fls), emit: ref_bwa_index_fls
        //path(ref_fasta_idx_fl), emit: ref_fasta_idx_fl
        //path(ref_fasta_fasta_dict), emit: ref_fasta_dct
        // -------------------------------------------------------------------


    script:
        base_name=sample_cram.baseName
        """
        bamcollate2 collate=${params.bamcollate_collate} \
            level=${params.bamcollate_level} \
            inputformat="cram" \
            filename=${sample_cram} \
            O="${base_name}.collated.bam"
        """
}
