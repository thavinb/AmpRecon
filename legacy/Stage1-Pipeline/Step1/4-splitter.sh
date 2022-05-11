samtools split --threads 4 -v \
--output-fmt cram,no_ref=1 \
-f $1/step1-4/%!.cram \
$1/step1-3.bam
