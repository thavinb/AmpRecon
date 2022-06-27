process bam_to_cram {
    /*
    * split BAM by read group into CRAM.
    */
    errorStrategy 'finish'
    publishDir "${params.results_dir}/${run_id}", overwrite: true

    input:
        tuple val(run_id), path(adapters_bam_file)

    output:
        tuple val(run_id), path("*.cram")

    script:
        """
        samtools split \
            --threads 4 \
            -v \
            --output-fmt "cram",no_ref=1 \
            -f "%!.cram" \
            ${adapters_bam_file}
        """
}
