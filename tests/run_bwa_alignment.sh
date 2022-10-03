#!/bin/bash


set -e 
cd tests
nextflow run test_align_bam.nf -c test_align_bam.config -profile standard  
