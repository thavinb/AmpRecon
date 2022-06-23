params.bamadapterclip_verbose = 1
params.bamadapterclip_level = 0

process clip_adapters {
    /*
    * removes identified adapters from bam
    */
    container ''

    input:
        tuple val(sample_tag), path(prealigned_bam)

    output:
        tuple val("${sample_tag}"), path("${base_name}.clipped.bam")

    script:
        base_name=prealigned_bam.baseName
        """
        bamadapterclip \
            verbose=${params.bamadapterclip_verbose} \
            level=${params.bamadapterclip_level} \
            < "${prealigned_bam}" \
            > "${base_name}.clipped.bam"
        """
}
