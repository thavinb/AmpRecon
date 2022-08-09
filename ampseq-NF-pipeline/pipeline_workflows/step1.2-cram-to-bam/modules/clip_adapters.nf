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
        //gambiarra alert --- this output is used on step1.2 only
        tuple val("${sample_tag}"), path("${base_name}.clipped.bam"), emit: tuple

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
// --- GAMBIARRA ALERT ---
// this output is used on step1.2 only. I was rushing to get step 1.3 running
// step 1.3 required some small changes to bam_to_fastq and align_bam had to be made,
// and the input for bam_to_fastq in step 1.2 comes from this module
