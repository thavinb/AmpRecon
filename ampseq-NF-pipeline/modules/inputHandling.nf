def validate_parameters() {
    // --- SANITY CHECKS ------------------------------------------------------
    def errors = 0
    // check if manifest file exists and if it was set
    if (params.manifest) {
        manifest_file=file(params.manifest)
        if (!manifest_file.exists()) {
            log.error("The manifest file specified does not exist.")
            errors += 1
        }
    }
    if (params.manifest==null){
        log.error("A manifest must be specified.")
        errors +=1
    }


    // check if reference fasta file exist and if it was set
    if (params.reference_fasta){
        reference_fasta = file(params.reference_fasta)
        if (!reference_fasta.exists()){
          log.error("The manifest file specified does not exist.")
          errors += 1
        }
    }
    if (params.reference_fasta==null){
        log.error("A reference fasta must be specified.")
        errors +=1
    }

    // check if output dir exists, if not create the default
    if (params.results_dir){
       results_path = file(params.results_dir)
       if (!results_path.exists()){
         log.warn("${results_path} does not exists, the dir will be created")
         results_path.mkdir()
       }
    }

    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}


def load_manifest_ch(){

  manifest_ch = Channel.fromPath(params.manifest) | splitCsv(header:true) | map {
                row-> tuple(row.run_id,
                            file(row.bcl_dir_path),
                            row.lane,
                            row.study_name,
                            row.read_group,
                            row.library)
                    }
  return manifest_ch
}
