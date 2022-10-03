#!/bin/bash


set -e 

nextflow run tests/test_align_bam.nf -c tests/test_align_bam.config -profile standard  
