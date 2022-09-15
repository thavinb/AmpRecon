#!/bin/bash


set -e


nextflow run tests/test_bamadapter_find.nf -c tests/test_bamadapter_find.config -profile standard 
wait
nextflow run tests/test_samtools_bam2cram.nf -c tests/test_samtools_bam2cram.config -profile standard   
wait
nextflow run tests/test_collate_alignments.nf -c tests/test_collate_alignments.config -profile standard
wait
nextflow run tests/test_bamreset.nf -c tests/test_bamreset.config -profile standard
wait
nextflow run tests/test_clip_adapters.nf -c tests/test_clip_adapters.config -profile standard
wait
nextflow run tests/test_bam_to_fastq.nf -c tests/test_bam_to_fastq.config -profile standard
wait
nextflow run tests/test_align_bam.nf -c tests/test_align_bam.config -profile standard
wait
nextflow run tests/test_scramble.nf -c tests/test_scramble.config -profile standard
