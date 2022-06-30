#!/bin/bash


set -e


nextflow run test_bamadapter_find.nf -c test_bamadapter_find.config -profile standard 
