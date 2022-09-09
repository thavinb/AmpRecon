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

  // a reference fasta is required for starting from 0 or 1.2a
  /*
  if (tag_provided=="0" || tag_provided=="1.2a") {
       check if reference fasta file exist and if it was set
      if (params.reference_fasta){
        reference_fasta = file(params.reference_fasta)
        if (!reference_fasta.exists()){
          log.error("The reference_fasta file specified ${params.reference_fasta} does not exist.")
          errors += 1
        }
      }
      if (params.reference_fasta==''){
        log.error("A reference fasta must be specified.")
        errors +=1
      }
  }
  */
  // --------------------------------------------------------------------------
  /*
  // redo_reference_fasta required checks
  if ((tag_provided=="0") || (tag_provided=="1.2a") ||(tag_provided="1.2b") || (tag_provided =="1.3")){
      // assert a new_reference is provided and exists
      if (params.redo_reference_fasta==null){
        log.error("A redo_reference fasta must be specified if start_from = ${tag_provided}.")
        errors +=1
      }
      if (params.redo_reference_fasta){
          redo_reference_fasta = file(params.redo_reference_fasta)
          if (!redo_reference_fasta.exists()){
            log.error("The redo_reference_fasta (${redo_reference_fasta}) file specified does not exist.")
            errors += 1
        }
      }
  }
  */
  if (params.irods_manifest && !(tag_provided== "1.2b")){
    log.warn("An iRODS manifest was provided (${params.irods_manifest}) but not needed for start_from = ${params.start_from}")
  }
  // --------------------------------------------------------------------------
  // A step1.1 out csv is required for step 1.2a
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
  // --------------------------------------------------------------------------
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
    // assert input_csv was provided and exists
    if (params.step1_3_in_csv==null){
      log.error("A 'step1_3_In_ch' must be specified if start_from = {$tag_provided}.")
      errors +=1
    }
    if (params.step1_3_in_csv){
      step1_3_in_csv = file(params.step1_3_in_csv)
      if (!step1_3_in_csv.exists()){
        log.error("The step1_3_in_csv file specified (${params.step1_3_in_csv}) does not exist.")
        errors += 1
      }
    }
  }

  // count errors and kill nextflow if any had been found
  if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
  }
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
          "1.2b":["1.2b","1.3"], // from iRODS
  ]
  // get the tag provided
  def tag_provided = params.start_from.toString()
  return steps_to_run_tagsMap.get(tag_provided)
}

