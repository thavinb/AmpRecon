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
        tuple val(tag), path(input_split_bam_file)
        path(input_merge_file)

    output:
        tuple val(tag), path("${output_file}")

    script:
        base_name=input_merge_file.baseName
        output_file="${base_name}.merged.bam"

        """
        bam12auxmerge \
            level=${params.bam12auxmerge_level} \
            rankstrip=${params.bam12auxmerge_rankstrip} \
            ranksplit=${params.bam12auxmerge_ranksplit} \
            zztoname=${params.bam12auxmerge_zztoname} \
            clipreinsert=${params.bam12auxmerge_clipreinsert} \
            "${input_merge_file}" \
            < "${input_split_bam_file}" \
            > "${output_file}"

        """
}
