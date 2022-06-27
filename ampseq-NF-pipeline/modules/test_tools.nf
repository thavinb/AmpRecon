process test_size {
    input:
        file(reference_file)
        file(test_file)

    output:
        stdout emit: exit_status

    script:
        """
        # Function to compare files
        function file_size {
            local test=\${1}
            local reference=\${2}

            absolute_test=\$(readlink -f \${test})
            absolute_ref=\$(readlink -f \${reference})
            size_test=\$(du \${absolute_test} | cut -f1)
            size_ref=\$(du \${absolute_ref} | cut -f1)

            if [[ "\${size_test}" == "\${size_ref}" ]]
            then
                echo "Passed. Size expected."
            else
                echo "Failed: Size not expected."
                exit 1
            fi
        }

        # Compare file size
        file_size ${test_file} ${reference_file}
        """
}

process test_binary {
    input:
        file(reference_file)
        file(test_file)

    output:
        stdout emit: exit_status

    script:
        """
        # Function to compare files
        function files_binary_equal {
            local test=\${1}
            local reference=\${2}

            local plat=\$(uname)
            if [[ "\${plat}" = "Darwin" ]]
            then
                absolute_test=\$( cd "\$(dirname "\${test}")" ; pwd -P )
                absolute_ref=\$( cd "\$(dirname "\${reference}")" ; pwd -P )
            else
                absolute_test=\$(readlink -f \${test})
                absolute_ref=\$(readlink -f \${reference})
            fi

            cmp \${absolute_test} \${absolute_ref}
            if [[ \$? == 0 ]]
            then
                echo "Passed. Files are binary equal."
            else
                echo "Failed: Files are not binary equal."
                exit 1
            fi
        }

        # Compare files for binary equality
        files_binary_equal ${test_file} ${reference_file}
        """
}

process test_value_equal {
    input:
        val value_1
        val value_2

    exec:
        assert value_1 == value_2
}


process subset_bam {

	//subsets 10% of reads from bam file for diff comparison

	input:
	tuple val(sample_id), path(bam_file)

	output:
	path("${output_file}")

	script:
	output_file="${sample_id}.test_subset.bam" 
	"""
	samtools view -bo ${output_file} -s 123.1 ${bam_file} 
	"""

}

process compare_bam_subset {

	//compares subset bams for equality 

	input:
	path(reference_bam)
	path(under_test)

	script:
	"""
	diff <(samtools view reference_bam) <(samtools view under_test) > differences.txt
	wait
	if [[ `wc -l < differences.txt` -ge 1 ]]; then
		echo subset bams differ 
		echo test failed
		exit 1
	fi
	"""

}


process check_md5sum {
    /**
    * Compares md5sum of input file with a given reference
    */

    input:
        file(input_file)
        val expected // reference md5sum

    script:
        """
        test_value=\$(md5sum ${input_file} | awk '{print \$1}')
        echo ${input_file} \${test_value} >&2
        if [[ "\${test_value}" != "${expected}" ]]; then
            echo "md5sum doesn't match" >&2
            exit 1
        fi
        """
}

process check_cram_md5sum {
    /**
    * Compares md5sum of input cram file with a given reference
    */

    input:
        path(input_file)
        val expected // reference md5sum

    script:
        """
        test_value=\$(samtools view -h ${input_file} | md5sum | awk '{print \$1}')
        echo ${input_file} \${test_value} >&2
        if [[ "\${test_value}" != "${expected}" ]]; then
            echo "md5sum doesn't match" >&2
            exit 1
        fi
        """
}

process check_vcf_md5sum {
    /**
    * Compares md5sum of input vcf file with a given reference
    */

    input:
        path(input_file)
        val expected // reference md5sum

    script:
        """
        test_value=\$(bcftools view -H ${input_file} | md5sum | awk '{print \$1}')
        echo ${input_file} \${test_value} >&2
        if [[ "\${test_value}" != "${expected}" ]]; then
            echo "md5sum doesn't match" >&2
            exit 1
        fi
        if [ \$(wc -l ${input_file}) -eq 0 ]; then
            echo "${input_file} is empty" >&2
            exit 1
        fi
        """
}


/**
* The two processes below take a tuple of file.baseName and two files as input
* It does this because when comparing two files, they need to be sorted by baseName
* to ensure that the correct files are being compared, the optimal way of sorting the files
* is to use the groupTuple() method, therefore the input must be a tuple.
*/

// check md5sum of two cram files
process compare_two_cram_files {
    input:
        tuple val(key), path(crams)

    output:
        stdout emit: exit_status

    script:
    """
        reference_value=\$(samtools view ${crams[0]} | md5sum | awk '{print \$1}')
        test_value=\$(samtools view ${crams[1]} | md5sum | awk '{print \$1}')
        if [[ "\${test_value}" != "\${reference_value}" ]]; then
            echo "md5sum doesn't match"
            exit 1
        fi

    """
}

// compare two vcf files
process compare_two_vcf_files {
        input:
        tuple val(key), path(vcfs)

        output:
        stdout emit: error_code

        script:
        """
        bcftools view --no-header ${vcfs[0]} > ref.vcf
        bcftools view --no-header ${vcfs[1]} > test.vcf
        diff ref.vcf test.vcf > vcf-diff.out
        rm vcf-diff.out
        """
}

