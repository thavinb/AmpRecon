bcftools filter \
--mode + \
--soft-filter LowDepth \
--exclude FORMAT/DP\<8 \
--output-type v \
< $1/step3-2/$2_1#$3.bcf \
> $1/step3-3/$2_1#$3.vcf
