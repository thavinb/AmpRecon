
def validate_parameters() {
  def errors = 0
  def valid_start_from_tags = ["0","1.2"]
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
  // a reference fasta is required for starting from 0 or 1.2
  if (tag_provided=="0" || tag_provided=="1.2") {
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

  // A step1.1 out csv is required for step 1.2
  if (params.start_from="1.2"){
      if (params.step1_2_in_csv==null){
        log.error("A step1_2_out_csv must be specified.")
        errors +=1
      }
      if (params.step1_2_in_csv){
        step1_2_in_csv = file(params.step1_2_in_csv)
        if (!step1_2_in_csv.exists()){
          log.error("The step1_2_in_csv ${params.step1_2_in_csv}file specified does not exist.")
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
