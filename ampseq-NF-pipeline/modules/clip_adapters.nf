params.bamadapterclip_verbose = 1
params.bamadapterclip_level = 0

process clip_adapters {
    /*
    * removes identified adapters from bam
    */

    input:
        val(sample_tag)
        path(prealigned_bam)

    output:
        val("${sample_tag}"), emit: sample_tag
        path("${base_name}.clipped.bam"), emit: clipped_bam

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
