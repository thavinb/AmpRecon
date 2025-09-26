// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process VALIDATE_MANIFEST {

    tag "${manifest.getBaseName()}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' : 
        'biocontainers/pandas:2.2.1' }"

    input:
        path(manifest)
        path(panel_setting)
        val(execution_mode)

    output:
        path("versions.yml"), emit: versions

    script:
        """
        validate_samplesheet.py ${manifest} ${panel_setting} ${execution_mode}

        cat <<-EOF > versions.yml
            ${task.process}:
                python: \$(python --version | cut -f2 -d ' ') 
        EOF
        """
}
