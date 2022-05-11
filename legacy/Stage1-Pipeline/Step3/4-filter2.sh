bcftools filter \
--mode + \
--soft-filter LowQual \
--exclude '%QUAL<15 || MQ<20' \
--output-type v \
< $1/step3-3/$2_1#$3.vcf \
> $1/step3-4/$2_1#$3.vcf
