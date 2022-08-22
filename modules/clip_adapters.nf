params.bamadapterclip_verbose = 1
params.bamadapterclip_level = 0

process CLIP_ADAPTERS {
    /*
    * removes identified adapters from bam
    */

    input:
        tuple val(run_id), path(prealigned_bam)

    output:
        //gambiarra alert --- this output is used on step1.2 only
        tuple val(run_id), path("${base_name}.clipped.bam")

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
