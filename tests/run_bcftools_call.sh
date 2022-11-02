#!/bin/bash


set -e 
cd tests
nextflow run test_bcftools_call.nf  -c test_bcftools.config -profile standard  
