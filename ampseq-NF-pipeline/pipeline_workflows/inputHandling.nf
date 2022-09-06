
def validate_parameters() {
  def errors = 0
  def valid_start_from_tags = ["0","1.2a","1.2b","1.3","S3"]
  // be sure that a start from is interpret as a string
  def tag_provided = params.start_from.toString()

  // --- SANITY CHECKS ------------------------------------------------------
  // check if output dir exists, if not create the default
  if (params.run_id == null){
      log.error("A run_id was not provided")
      errors += 1
  }

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
      // check if all expected input parameters was set and if the file exists
    if (params.study_name == null){
        log.error("A run_id was not provided.")
        errors += 1
    }

    if (params.lane == null) {
      log.error("A value for parameter lane was not provided.")
      errors +=1
    }

    if (params.read_group == null){
      log.error("A value for parameter read_group was not provided.")
      errors +=1
    }

    if (params.library == null){
      log.error("A value for parameter library was not provided.")
      errors +=1
    }

    if (params.bcl_dir == null) {
      log.error("A value for parameter library was not provided.")
      errors +=1
    }
    if (!(params.bcl_dir == null)) {
      bcl_dir_path = file(params.bcl_dir)
      if (!bcl_dir_path.exists()){
        log.error("${bcl_dir_path} does not exists.")
        errors += 1
      }
    }
  }
  // the input csv is only needed if starts from 0, otherwise should be ignored
  if (!(tag_provided== "0") && !(tag_provided=="S3")){
    if (params.study_name){
     log.warn("A study_name was provided (${params.study_name}) but ignored (not needed for start_from = ${params.start_from})")
    }

    if (params.lane){
     log.warn("A lane was provided (${params.lane}) but ignored (not needed for start_from = ${params.start_from})")
    }
    if (params.read_group){
     log.warn("A read_group was provided (${params.read_group}) but ignored (not needed for start_from = ${params.start_from})")
    }
    if (params.library){
     log.warn("A library was provided (${params.library}) but ignored (not needed for start_from = ${params.start_from})")
    }
    if (params.bcl_dir){
     log.warn("A bcl_dir was provided (${params.bcl_dir}) but ignored (not needed for start_from = ${params.start_from})")
    }
  }

  // --------------------------------------------------------------------------

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
          "S3":["S3","1.2a", "1.3"], // in country from S3
  ]
  // get the tag provided
  def tag_provided = params.start_from.toString()
  return steps_to_run_tagsMap.get(tag_provided)
}
