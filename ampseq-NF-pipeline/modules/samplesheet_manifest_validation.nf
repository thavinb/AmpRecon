process validate_samplesheet_manifest {

    input:
        tuple val(run_id), path(manifest)
    output:
        tuple val(run_id), path("${manifest}")

    script:
    """
    python ${workflow.projectDir}/modules/manifest_validation.py -m ${manifest}
    """
    }
