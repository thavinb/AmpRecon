
process mpileup {

	input:
	tuple val(sample_tag), file(input_bam)
	file(reference_fa)
	file(annotation_vcf)


	output:
	tuple val(sample_tag), path("*.vcf")


	script:
	"""
	bcftools mpileup \
	--min-BQ 20 \
	--annotate FORMAT/AD,FORMAT/DP \
	--max-depth 50000 \
	--targets-file ${annotation_vcf} \
	--fasta-ref ${reference_fa} \
	--output-type u \
	${input_bam} > ${sample_tag}.vcf
	"""

}

