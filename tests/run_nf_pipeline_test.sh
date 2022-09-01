#!/bin/bash


set -e

nextflow run tests/test_pipeline_in-country_setup.nf -c tests/test_pipeline_setup.config -profile standard 
wait
nextflow run ampseq-NF-pipeline/main.nf -c tests/test_pipeline_in-country_entry_point.config -profile sanger_default
wait
nextflow run tests/validate_pipeline_in-country_entry_point.nf -c tests/test_pipeline_in-country_entry_point.config -profile sanger_default
#wait
#nextflow run ampseq-NF-pipeline/main.nf -c tests/test_pipeline_iRODS_entry_point.config -profile sanger_default
#wait
#nextflow run tests/validate_pipeline_iRODS_entry_point.nf -c tests/test_pipeline_iRODS_entry_point.config -profile sanger_default
