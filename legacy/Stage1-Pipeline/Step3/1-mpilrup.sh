bcftools mpileup \
--min-BQ 20 \
--annotate FORMAT/AD,FORMAT/DP \
--max-depth 50000 \
--targets-file $4.annotation.vcf \
--fasta-ref $4.fasta \
--output-type u \
$1/step2-12/$2_1#$3.bam \
> $1/step3-1/$2_1#$3.bcf
