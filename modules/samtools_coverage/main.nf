// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process SAMTOOLS_COVERAGE {
    /*
    * Get coverage from bam file
    */
    tag "${meta.uuid}"
    label 'samtools'
    errorStrategy 'ignore'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.coverage"), emit: coverage
        path "versions.yml"                , emit: versions

    script:
        def prefix   = "${meta.uuid}"
        """
        if samtools coverage --help > /dev/null 2>&1; then
            samtools coverage ${bam} > ${prefix}.coverage
        else
            samtools depth -a ${bam} | awk '
                BEGIN {
                    print "Chrom\\tLength\\tCovbases\\tCoverage(%)\\tMeanDepth"
                }
                {
                    if (\$1 != prev && NR > 1) {
                        printf "%s\\t%d\\t%d\\t%.2f\\t%.2f\\n", prev, len, covered, (covered/len)*100, sum/len
                        sum=0; covered=0; len=0
                    }

                    prev = \$1
                    len++
                    sum += \$3
                    if (\$3 > 0) covered++
                }
                END {
                    if (len > 0) {
                         printf "%s\\t%d\\t%d\\t%.2f\\t%.2f\\n", prev, len, covered, (covered/len)*100, sum/len
                    }
                }' > ${prefix}.coverage
        fi

        cat <<-EOF > versions.yml
        "${task.process}":
            samtools: "\$( samtools --version | grep -E "^samtools" | cut -f2 -d ' '  )"
            htslib: "\$( samtools --version | grep -E '^Using htslib' | cut -f 3 -d ' ')"
        EOF
        """
}
