bcftools call \
--multiallelic-caller \
--keep-alts \
--skip-variants indels \
--ploidy-file $4.ploidy \
--output-type u \
< $1/step3-1/$2_1#$3.bcf \
> $1/step3-2/$2_1#$3.bcf
