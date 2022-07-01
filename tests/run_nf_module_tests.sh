#!/bin/bash


set -e


nextflow run test_bamadapter_find.nf -c test_bamadapter_find.config -profile standard 
wait
nextflow run test_samtools_bam2cram.nf -c test_samtools_bam2cram.config -profile standard   
