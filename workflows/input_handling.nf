def validate_parameters() {

  def errors = 0
  def valid_execution_modes = ["in-country", "irods"]
  
  // --- SANITY CHECKS ------------------------------------------------------
  // check if output dir exists, if not create the default
  if (params.results_dir){
       results_path = file(params.results_dir)
       if (!results_path.exists()){
         log.warn("${results_path} does not exists, the dir will be created")
         results_path.mkdir()
       }
  }

  // check if execution mode is valid
  if (!valid_execution_modes.contains(params.execution_mode)) {
    log.error("The execution mode provided (${params.execution_mode}) is not valid. valid modes = ${valid_execution_modes}")
    errors += 1
  }

  // check if all params required for in country were provided
  if (params.execution_mode == "in-country"){
    if (params.run_id == null){
      log.error("A run_id parameter must be provided for execution mode '${params.execution_mode}'.")
      errors += 1
    }

    if (params.bcl_dir == null && params.s3_launch_uuid == null){
      log.error("Either a bcl directory or a uuid corresponding to a run directory in s3 must be specified.")
      errors += 1
    } 
    

    if (params.lane == null){
      log.error("A lane parameter must be provided for execution mode '${params.execution_mode}'.")
      errors += 1
    }
    
    if (params.study_name == null){
      log.error("A study_name parameter must be provided for execution mode '${params.execution_mode}'.")
      errors += 1
    }
    if (params.read_group == null){
      log.error("A read_group parameter must be provided for execution mode '${params.execution_mode}'.")
      errors += 1
    }
    
    if (params.library == null){
      log.error("A library parameter must be provided for execution mode '${params.execution_mode}'.")
      errors += 1
    }

  }

  // check if all params required for irods were provided
  if (params.execution_mode == "irods"){
    if (params.irods_manifest == null){
      log.error("An irods_manifest parameter must be provided for execution mode '${params.execution_mode}'.")
      errors += 1
    }
    if (params.irods_manifest){
      irods_manifest = file(params.irods_manifest)
      if (!irods_manifest.exists()){
        log.error("The irods manifest file specified (${params.irods_manifest}) does not exist.")
        errors += 1
      }
    }
  }

  // check if all s3 required parameters were provided
  if (!(params.s3_launch_uuid==null)){
    if (params.s3_bucket_input == null){
      log.error("A s3_bucket_input parameter must be provided if download_from_s3 is set to '${params.download_from_s3}'.")
      errors += 1
    } 
  }

  if (params.upload_to_s3){
    if (params.s3_bucket_output == null){
      log.error("A s3_bucket_output parameter must be provided if download_from_s3 is set to '${params.upload_from_s3}'.")
      errors += 1
    }
  }
  // count errors and kill nextflow if any had been found
  if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
  }
}
