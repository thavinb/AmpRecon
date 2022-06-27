#!/bin/bash


#run bambi i2b outside nextflow environment to generate test files to validate against


#BAMBI I2B

singularity exec --bind /lustre,/nfs,/software /lustre/scratch123/pam/dt5/personal/malpipe2/ss52-ampseq/ampseq-pipeline/containers/bambi.sif bambi i2b --intensity-dir=./21045/Data/Intensities --basecalls-dir=./21045/Data/Intensities/BaseCalls --lane="1" --read-group-id="21045_1" --study-name="Test" --threads=8 --output-file=21045_bambi_i2b.bam --output-fmt="bam" --compression-level=9   
