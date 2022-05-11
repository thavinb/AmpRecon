python3 7-mergeHeaders.py $3 $4 $1 $2 | samtools reheader - \
$1/step2-6/$2_1#$3.bam \
> $1/step2-7/$2_1#$3.bam
