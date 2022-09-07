params.bamcollate_collate = 1
params.bamcollate_level = 9

process collate_alignments {
    /*
    * collates bam/cram reads/alignments by name
    */
    errorStrategy 'ignore'

    //publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        //tuple val(run_id), path(cram_fl), val(sample_tag)
        val(run_id)
        path(sample_cram)
        val(sample_tag)

    output:
        val("${sample_tag}"), emit: sample_tag
        path("${base_name}.collated.bam"), emit: collated_bam

    script:
        base_name=sample_cram.simpleName
        """
        bamcollate2 collate=${params.bamcollate_collate} \
            level=${params.bamcollate_level} \
            inputformat="cram" \
            filename=${sample_cram} \
            O="${base_name}.collated.bam"
        """
}
