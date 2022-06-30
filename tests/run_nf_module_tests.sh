#!/bin/bash


set -e


nextflow run ./test_bamadapter_find.nf -c ../ampseq-NF-pipeline/nextflow.config 
