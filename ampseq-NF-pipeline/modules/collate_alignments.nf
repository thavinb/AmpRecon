params.bamcollate_collate = 1
params.bamcollate_level = 9

process collate_alignments {
    /*
    * collates bam/cram reads/alignments by name
    */
    container ''

    input:
        val(tag)
        path(input_file)

    output:
        tuple val("${tag}"), path("${base_name}.collated.bam")

    script:
        base_name=input_file.baseName
        """
        bamcollate2 collate=${params.bamcollate_collate} \
            level=${params.bamcollate_level} \
            inputformat="cram" \
            filename=${input_file} \
            O="${base_name}.collated.bam"
        """
}
