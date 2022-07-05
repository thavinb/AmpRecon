#!/bin/bash


singularity exec --bind /lustre,/nfs,/software /lustre/scratch123/pam/dt5/personal/malpipe2/ss52-ampseq/ampseq-pipeline/containers/samtools_1.15.sif samtools split --threads 4 -v --output-fmt "cram",no_ref=1 -f "%!.cram" 21045.adapters.bam 
