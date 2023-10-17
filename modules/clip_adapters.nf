// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

params.bamadapterclip_verbose = 1
params.bamadapterclip_level = 0

process clip_adapters {
    /*
    * Removes identified adapters from bam
    */

    input:
        tuple val(file_id), path(input_bam)

    output:
        tuple val("${file_id}"), path("${base_name}.clipped.bam")

    script:
        base_name=input_bam.simpleName
        """
        bamadapterclip \
            verbose=${params.bamadapterclip_verbose} \
            level=${params.bamadapterclip_level} \
            < "${input_bam}" \
            > "${base_name}.clipped.bam"
        """
}
