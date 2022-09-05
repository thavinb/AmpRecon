
def validate_parameters() {
  def errors = 0
  def valid_start_from_tags = ["0","1.2a","1.2b","1.3"]
  // be sure that a start from is interpret as a string
  def tag_provided = params.start_from.toString()

  // --- SANITY CHECKS ------------------------------------------------------
  // check if output dir exists, if not create the default
  if (params.results_dir){
       results_path = file(params.results_dir)
       if (!results_path.exists()){
         log.warn("${results_path} does not exists, the dir will be created")
         results_path.mkdir()
       }
  }
  if (!valid_start_from_tags.contains(tag_provided)) {
    log.error("The start_from provided (${params.start_from}) is not valid.")
    errors += 1
  }
  // ------------------------------------------------------------------------
  // Check if everything required was provided for a given entry point
  if (tag_provided == "0"){
      // check if manifest was set and if the file exists
      if (params.manifest) {
        manifest_file=file(params.manifest)
        if (!manifest_file.exists()) {
            log.error("The manifest file specified does not exist.")
            errors += 1
          }
      }
      if (params.manifest==''){
        log.error("A manifest must be specified.")
        errors +=1
      }

  }
  // the input csv is only needed if starts from 0, otherwise should be ignored
  if (params.manifest && !(tag_provided== "0")){
    log.warn("A manifest was provided (${params.manifest}) but ignored (not needed for start_from = ${params.start_from})")
  }

  // a reference fasta is required for starting from 0 or 1.2
  if (tag_provided=="0" || tag_provided=="1.2a") {
      // check if reference fasta file exist and if it was set
      if (params.reference_fasta){
        reference_fasta = file(params.reference_fasta)
        if (!reference_fasta.exists()){
          log.error("The manifest file specified does not exist.")
          errors += 1
        }
      }
      if (params.reference_fasta==''){
        log.error("A reference fasta must be specified.")
        errors +=1
      }
  }
  if (params.irods_manifest && !(tag_provided== "1.2b")){
    log.warn("An iRODS manifest was provided (${params.irods_manifest}) but not needed for start_from = ${params.start_from}")
  }

  // A step1.1 out csv is required for step 1.2
  if (tag_provided=="1.2a"){
      if (params.step1_2_in_csv==''){
        log.error("A step1_2_out_csv must be specified.")
        errors +=1
      }
      if (params.step1_2_in_csv){
        step1_2_in_csv = file(params.step1_2_in_csv)
        if (!step1_2_in_csv.exists()){
          log.error("The step1_2_in_csv (${params.step1_2_in_csv}) file specified does not exist.")
          errors += 1
        }
      }
      // the input csv is only needed if starts from 0, otherwise should be ignored
  }

  // an irods manifest for tag 1.2b (iRODS pulling BAM files)
  if (tag_provided=="1.2b"){
    if (params.irods_manifest==""){
      log.error("A irods_manifest must be specified.")
      errors +=1
    }
    if (params.irods_manifest){
      irods_manifest = file(params.irods_manifest)
      if (!irods_manifest.exists()){
        log.error("The irods manifest file specified (${params.irods_manifest}) does not exist.")
        errors += 1
      }
    }
  }
  // --- 1.3 tag checks -------------------------------------------------------
  if (tag_provided=="1.3"){

  }

  // count errors and kill nextflow if any had been found
  if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
  }
}


def load_manifest_ch(){

  manifest_ch = Channel.fromPath(params.manifest) | splitCsv(header:true) |
                map {row-> tuple(row.run_id,
                                 row.bcl_dir_path,
                                 row.lane,
                                 row.study_name,
                                 row.read_group,
                                 row.library)
                    }
  return manifest_ch
}

def load_steps_to_run(){
  /*
  this function define which tags / subworkflows should be considered for
  current run.
  INPUT: params.start_from
  */
  // the steps_to_run_tags list should be ordered according to expected
  // execution order
  // PS: for now this is okay, but if we need subworkflows that create branchs
  //     to a set of exclusive process, we may have to change the current structure
  //     behaviour

  // define a set tags to steps to run map
  def steps_to_run_tagsMap = [
          "0":["0","1.2a", "1.3"], // in country by default
          "1.2a":["1.2a","1.3"],
          "1.2b":["1.2b","1.3"] // from iRODS
  ]
  // get the tag provided
  def tag_provided = params.start_from.toString()
  return steps_to_run_tagsMap.get(tag_provided)
}
