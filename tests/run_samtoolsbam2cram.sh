#!/bin/bash


set -e

nextflow run tests/test_samtools_bam2cram.nf -c tests/test.config -profile standard   
