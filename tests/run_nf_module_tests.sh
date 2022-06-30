#!/bin/bash


set -e


nextflow run tests/test_bamadapter_find.nf -c tests/test_bamadapter_find.config -profile standard 
