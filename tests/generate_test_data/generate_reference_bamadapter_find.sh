#!/bin/bash


singularity exec --bind /lustre,/nfs,/software /lustre/scratch123/pam/dt5/personal/malpipe2/ss52-ampseq/ampseq-pipeline/containers/biobambam2.sif bamadapterfind level=9 < 21045.test_i2b.subset.bam > 21045.adapters.bam 
