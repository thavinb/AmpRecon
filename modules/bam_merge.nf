params.bam12auxmerge_level = 0
params.bam12auxmerge_rankstrip = 1
params.bam12auxmerge_ranksplit = 0
params.bam12auxmerge_zztoname = 0
params.bam12auxmerge_clipreinsert = 1


process bam_merge {
    /*
    * Merges BAM files containing same set of reads
    */
    container ''

    input:
        tuple val(sample_tag), path(reheadered_bam), path(split_bam)

    output:
        val(sample_tag), emit: sample_tag
        path("${merged_bam}"), emit: merged_bam

    script:
        base_name=reheadered_bam.simpleName
        merged_bam="${base_name}.merged.bam"

        """
        bam12auxmerge \
            level=${params.bam12auxmerge_level} \
            rankstrip=${params.bam12auxmerge_rankstrip} \
            ranksplit=${params.bam12auxmerge_ranksplit} \
            zztoname=${params.bam12auxmerge_zztoname} \
            clipreinsert=${params.bam12auxmerge_clipreinsert} \
            "${reheadered_bam}" \
            < "${split_bam}" \
            > "${merged_bam}"
        """
}
