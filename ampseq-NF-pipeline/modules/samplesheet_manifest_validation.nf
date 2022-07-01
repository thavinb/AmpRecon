process validate_samplesheet_manifest {

    input:
        tuple val(run_id), path(manifest)
    output:
        tuple val(run_id), path("${manifest}")

    script:
    """
    python3 ${workflow.projectDir}/modules/manifest_validation.py -m ${manifest}
    """
    }
