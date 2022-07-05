#!/bin/bash



singularity exec --bind /lustre,/nfs,/software /lustre/scratch123/pam/dt5/personal/malpipe2/ss52-ampseq/ampseq-pipeline/containers/bambi.sif bambi decode --input-fmt="bam" --metrics-file=21045.subset.metrics --barcode-file=../tags.tsv --output=21045_bambi_decode.subset.bam  --output-fmt="bam" --compression-level=0 21045.test_i2b.subset.bam     
