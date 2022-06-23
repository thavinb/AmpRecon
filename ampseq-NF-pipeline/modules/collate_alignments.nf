params.bamcollate_collate = 1
params.bamcollate_level = 9

process collate_alignments {
    /*
    * collates bam/cram reads/alignments by name
    */
    errorStrategy 'ignore'

    publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

    input:
        tuple val(run_id), path(sample_cram), val(sample_tag)
        //val(tag) // be sure what tag supose to mean here (not run id a supose)
        //path(input_file)

    output:
        tuple val("${sample_tag}"), path("${base_name}.collated.bam")

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
