/*
    | DESIGNATE_PANEL_RESOURCES |-----------------------------------------
    
    This workflow assigns each sample_tag the appropriate set of
    reference files, reference indexes and auxiliary files

    It keeps the manipulation of the sample tags and reference 
    files channel all in one place and forces the input to be a 
    specific size, simplifying the addition of new columns
    ------------------------------------------------------------------
*/

workflow DESIGNATE_PANEL_RESOURCES {
    take:
        sample_tag_ch // tuple (new_sample_id, panel_name)
        reference_ch // tuple (fasta, panel_name, snp_list)
    
    main:
        sample_tag_ch
          |  combine(reference_ch,  by: 1) // tuple (panel_name, new_sample_id, fasta,snp_list)
          |  map{it -> tuple(it[1], it[0], it[2], it[3])}
          |  set{sample_tag_reference_files_ch}
        // tuple (new_sample_id, panel_name, path/to/reference/genome, snp_list)

    emit:
        sample_tag_reference_files_ch // tuple (new_sample_id, panel_name, path/to/reference/genome, snp_list)
}
