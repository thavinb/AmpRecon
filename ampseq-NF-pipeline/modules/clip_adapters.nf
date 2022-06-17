params.bamadapterclip_verbose = 1
params.bamadapterclip_level = 0

process clip_adapters {
    /*
    * removes identified adapters from bam
    */
    container ''

    input:
        tuple val(tag), path(input_file)

    output:
        tuple val("${tag}"), path("${base_name}.clipped.bam")

    script:
        base_name=input_file.baseName
        """
        bamadapterclip \
            verbose=${params.bamadapterclip_verbose} \
            level=${params.bamadapterclip_level} \
            < "${input_file}" \
            > "${base_name}.clipped.bam"
        """
}
