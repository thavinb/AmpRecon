bambi select \
--compression-level=0 \
--input $1/step2-1/$2_1#$3.cram,$1/step2-9/$2_1#$3.bam \
--output $1/step2-10/$2_1#$3-phix.bam,$1/step2-10/$2_1#$3.bam \
-m $1/step2-10/$2_1#$3_bam_alignment_filter_metrics.json
