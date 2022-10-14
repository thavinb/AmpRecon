#!/bin/bash


set -e
cd tests
nextflow run test_samtools_bam2cram.nf -c test.config -profile standard   
