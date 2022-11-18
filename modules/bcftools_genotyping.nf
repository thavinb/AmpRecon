process bcftools_mpileup {
    /*
    * Creates an uncompressed BCF file containing calculated genotype likelihoods for every possible genomic position supported by the BAM
    */
    label 'bcftools'
    input:
        tuple val(sample_tag), path(input_bam), path(input_bam_index), val(reference_file), val(reference_annotation_vcf)

    output:
        tuple val(sample_tag), path("${output_bcf}")

    script:
        base_name = input_bam.baseName
        output_bcf="${base_name}.bcf"

        """
	bcftools mpileup \
	--min-BQ 20 \
	--annotate FORMAT/AD,FORMAT/DP \
	--max-depth 50000 \
	--targets-file "${reference_annotation_vcf}" \
	--fasta-ref "${reference_file}" \
	--output-type u \
	"${input_bam}" \
	> "${output_bcf}"
        """
}

process bcftools_call {
    /*
    * Calls SNPs from a BCF file containing all possible genotype likelihoods across genome
    */
    label 'bcftools'
    input:
        tuple val(sample_tag), path(input_bcf)

    output:
        tuple val(sample_tag), path("${output_bcf}")

    script:
        base_name = input_bcf.baseName
        output_bcf="${base_name}.bcftools_genotyped.bcf"
        ploidy="* * * * 2"

        """
        echo "${ploidy}" > ploidy_file.ploidy

	bcftools call \
	--multiallelic-caller \
	--keep-alts \
	--skip-variants indels \
	--ploidy-file "ploidy_file.ploidy" \
	--output-type u \
	< "${input_bcf}" \
	> "${output_bcf}"
        """
}

process bcftools_filter {
    /*
    * SNPs in the input BCF file are filtered and output as an uncompressed VCF file
    */
    publishDir "${params.results_dir}/", overwrite: true
    label 'bcftools'
    input:
        tuple val(sample_tag), path(input_bcf)

    output:
        tuple val(sample_tag), path("${output_vcf}")

    script:
        base_name = input_bcf.baseName
        intermediate_vcf="${base_name}.intermediate.vcf"
        output_vcf="${base_name}.vcf"

	// Had to escape backslash character in FORMAT line of script

        """
	bcftools filter \
	--mode + \
	--soft-filter LowDepth \
	--exclude FORMAT/DP\\<8 \
	--output-type v \
	< "${input_bcf}" \
	> "${intermediate_vcf}"

	bcftools filter \
	--mode + \
	--soft-filter LowQual \
	--exclude '%QUAL<15 || MQ<20' \
	--output-type v \
	< "${intermediate_vcf}" \
	> "${output_vcf}"

        rm "${intermediate_vcf}"
        """
}

