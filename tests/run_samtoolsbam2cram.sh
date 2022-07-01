#!/bin/bash


set -e

nextflow run tests/test_samtools_bam2cram.nf -c tests/test_samtools_bam2cram.config -profile standard   
