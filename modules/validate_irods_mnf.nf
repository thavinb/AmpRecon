#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process validate_irods_mnf {
    label 'pythonBox'
    input:
        path(irods_mnf_file)
        path(pannel_settings_file)

    script:
    """
    validate_irods_mnf.py ${irods_mnf_file} ${pannel_settings_file}
    """
}