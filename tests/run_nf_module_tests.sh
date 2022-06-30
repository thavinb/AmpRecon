#!/bin/bash


set -e


nextflow run tests/test_bamadapter_find.nf -c ampseq-NF-pipeline/nextflow.config -profile standard 
