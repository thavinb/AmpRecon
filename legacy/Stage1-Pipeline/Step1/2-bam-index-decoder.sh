bambi decode --output=$1/step1-2.bam \
 --compression-level=0 \
 --metrics-file=$2/metrics.txt \
 --barcode-file=$3 \
 $1/step1-1.bam
