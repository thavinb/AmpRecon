//params.pannels_settings = "/Users/ad45/Projects/git_stuff/ampseq-pipeline/ampseq-NF-pipeline/sandbox/pannels_settings.csv"

def validatePannelSettings(row){
    
    valid_pannel_names = ["PFA_GRC1_v1.0","PFA_GRC2_v1.0","PFA_Spec"]
    def errors = 0 
    // check if pannel name is valid
    if (!valid_pannel_names.contains(row.pannel_name)){
        log.error("${row.pannel_name} is not valid. Valid pannel names:\n ${valid_pannel_names}")
        errors += 1
    }

    // check if align_to columns is a valid path
    align_to_path = file(row.aligns_to)
    if (!align_to_path.exists()){
        log.error("${row.aligns_to} provided for ${row.pannel_name} does not exist.")
        errors += 1
    }

    // TODO: check if all expected files are present on the resource bundle provided
    
    // if any errors, kill the pipeline run
    
    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

workflow PARSE_PANNEL_SETTINGS {

    // if a csv was provided, mount channels for a given resource bundle
    if (!(params.pannels_settings==null)){
    
        Channel.fromPath(params.pannels_settings, checkIfExists: true)
            | splitCsv(header: true, sep: ',')
            | map { row ->
                        // validate row settings
                        validatePannelSettings(row)
                        // TODO: validate if expected files were provided
                        tuple(
                            file("${row.aligns_to}/*.fasta"),
                            row.pannel_name,
                            file("${row.aligns_to}/*.fasta.*")
                            )
                }
            | set { reference_ch }
    }
    // if pannels settings csv is not provided, just use the files on the repo
    else {
        reference_ch = Channel.from(
            [file("${params.reference_dir}/grc1/*.fasta"), "PFA_GRC1_v1.0" , file("${params.reference_dir}/grc1/*.fasta.*")],
            [file("${params.reference_dir}/grc2/*.fasta"), "PFA_GRC2_v1.0", file("${params.reference_dir}/grc2/*.fasta.*")],
            [file("${params.reference_dir}/spec/*.fasta"), "PFA_Spec", file("${params.reference_dir}/spec/*.fasta.*")]
        )
    }

    emit:
        reference_ch
}
