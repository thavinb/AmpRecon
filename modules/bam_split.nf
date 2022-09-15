params.bam12split_verbose = 0
params.bam12split_level = 0

process bam_split {
    /*
    * Splits BAM rank pairs to single ranks per read
    */
    label 'biobambam2'


    input:
        tuple val(sample_tag), path(reheadered_bam)

    output:
        tuple val(sample_tag), path("*.split.bam")

    script:
      base_name = reheadered_bam.getBaseName()
      """
      bam12split verbose=${params.bam12split_verbose} level=${params.bam12split_level} \
          < "${reheadered_bam}" \
          > "${base_name}.split.bam"
      """
}
