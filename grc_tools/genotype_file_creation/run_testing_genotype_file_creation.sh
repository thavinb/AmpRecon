#!/bin/bash


set -e 
cd tests
nextflow run test_genotype_file_creation.nf -c test_genotype_file_creation.config -profile standard
