// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process split_bam_by_readgroup {
    /*
    * split BAM by read group into CRAM.
    */

    input:
        path(adapters_bam_file)

    output:
        path("*.cram")

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

