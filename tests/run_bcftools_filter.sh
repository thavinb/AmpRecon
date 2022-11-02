#!/bin/bash


set -e 
cd tests
nextflow run test_bcftools_filter.nf  -c test_bcftools.config -profile standard  
