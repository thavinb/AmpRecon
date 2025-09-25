//
// Subworkflow contains useful utility functions and workflows
//

include { VALIDATE_MANIFEST } from '../../modules/validate_manifest/main'

/*
-----------------------------------------------------------------------------------
    SUBWORKFLOW FOR INITIALISE PIPELINE
-----------------------------------------------------------------------------------
*/

workflow PIPELINE_INIT {
    take: 
        help           // boolean: Display help and exit
        monochrome_logs   // boolean: Do not use coloured log outputs
        outdir            //  string: The output directory where the results will be saved
        input             //  string: Path to input samplesheet

    main:
        // Show help message
        if (help) {
           printHelp()
           exit 0
        }       
        // Print params before running
        loggingInit(monochrome_logs)

        // Validate All Input Params
        validateInputParams()

        // Validate and create referecne channel
        // tuple ( ref.fasta, snp )
        reference_ch = validateReferencePanels(params.panels_settings)

        if (params.execution_mode == "cram") {

            // validate sample sheet 
            VALIDATE_MANIFEST(
				params.manifest,
				params.panels_settings,
				params.execution_mode
			)

            // create irods channel 
            // tuple ( [uuid, id, panel, genome, snps], cram_path )
            irods_ch = Channel.fromPath("${params.manifest}")
                        | splitCsv(header: true, sep: '\t')
                        | map { row -> 
                            def meta = [
                                uuid: row.sample_id + "_" + row.primer_panel,
                                id: row.sample_id,
                                panel: row.primer_panel
                            ]
                            tuple( row.primer_panel, meta, row.cram_path ) }
                        | combine( reference_ch, by:0 )
                        | map { it ->  groupPanelReference(it) }
                        | set { input_ch }
        }

        if (params.execution_mode == "fastq") {

            // validate sample sheet
            VALIDATE_MANIFEST(
				params.manifest,
				params.panels_settings,
				params.execution_mode
			)

            // create fastq channel
            // tuple ( [uuid, id, panel, genome, snps], [fastq1_path, fastq2_path] )
            fastq_ch = Channel.fromPath("${params.manifest}")
                        | splitCsv(header: true, sep: '\t')
                        | map { row -> 
                            def meta = [
                                uuid: row.sample_id + "_" + row.primer_panel,
                                id: row.sample_id,
                                panel: row.primer_panel
                            ]
                            tuple( row.primer_panel, meta, [row.fastq1_path, row.fastq2_path] ) }
                        | combine( reference_ch, by:0 )
                        | map { it ->  groupPanelReference(it) }
                        | set { input_ch }
        }

    emit: 
        input_ch
}

/*
-----------------------------------------------------------------------------------
    SUBWORKFLOW FOR PIPELINE COMPLETION
-----------------------------------------------------------------------------------
*/
workflow PIPELINE_COMPLETION {
}
/*
-----------------------------------------------------------------------------------
    UTILS FUNCTION
-----------------------------------------------------------------------------------
*/

def loggingInit(monochrome_logs) {
    // logging info ----------------------------------------------------------------
    // This part of the code is based on the one present at FASTQC PIPELINE (https://github.com/angelovangel/nxf-fastqc/blob/master/main.nf)

    /*
    * ANSI escape codes to color output messages
    */
    if (!params.monochrome_logs) {
        ANSI_GREEN = "\033[1;32m"
        ANSI_RED   = "\033[1;31m"
        ANSI_RESET = "\033[0m"
    } else {
        ANSI_GREEN = ""
        ANSI_RED   = ""
        ANSI_RESET = ""
    }

    log.info"""
            ===========================================
             AMPRECON ${workflow.manifest.version}
             Used parameters:
            -------------------------------------------
             Pipeline:
             --execution_mode           : ${params.execution_mode}
             --panels_settings          : ${params.panels_settings}
             --containers_dir           : ${params.containers_dir}
             --batch_id                 : ${params.batch_id}

             Manifest:
             --manifest                 : ${params.manifest}

             Results Directory:
             --results_dir              : ${params.results_dir}

             GRC:
             --grc_settings_file_path   : ${params.grc_settings_file_path}
             --chrom_key_file_path      : ${params.chrom_key_file_path}
             --kelch_reference_file_path: ${params.kelch_reference_file_path}
             --codon_key_file_path      : ${params.codon_key_file_path}
             --drl_information_file_path: ${params.drl_information_file_path}
             --no_plasmepsin            : ${params.no_plasmepsin}
             --no_kelch                 : ${params.no_kelch}

             (DEBUG):
             --DEBUG_tile_limit         : ${params.DEBUG_tile_limit}
             --DEBUG_takes_n_bams       : ${params.DEBUG_takes_n_bams}
             --DEBUG_no_coi             : ${params.DEBUG_no_coi}
            ------------------------------------------
             Runtime data:
            -------------------------------------------
             Running with profile       : ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
             Running as user            : ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
             Launch dir                 : ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
             Base dir                   : ${ANSI_GREEN}${baseDir}${ANSI_RESET}
             ------------------------------------------
             """.stripIndent()
}

def printHelp() {
  log.info """
  AMPRECON ${workflow.manifest.version}
  Usage:
    (cram):
    nextflow run main.nf \\
      -profile standard \\
      -c conf/pfalciparum.config \\
      --execution_mode cram \\
      --batch_id 21045 \\
      --manifest manifest.tsv \\
      --containers_dir ./containers_dir/ \\
      --results_dir output

    (fastq):
    nextflow run main.nf \\
      -profile standard \\
      -c conf/pfalciparum.config \\
      --execution_mode cram \\
      --batch_id 21045 \\
      --manifest manifest.tsv \\
      --containers_dir ./containers_dir/ \\
      --results_dir output

  Description:
    AmpRecon is a bioinformatics analysis pipeline for amplicon sequencing data.
    Currently supporting alignment and SNP variant calling on paired-end Illumina sequencing data.

    *for a complete description of input files and parameters check the README file in the code repository

  Options:
    Inputs:
      (required)
      --execution_mode : sets the entry point for the pipeline ("cram" or "fastq")
      --manifest       : path to the manifest file
      --batch_id       : id to be used for the batch of data to be processed. 
                         This ID is only used to prefix output files and readgroup names in cram files.
                         It can be a run ID or any other identifier that makes sense for your data.

      (grc_creation)
      --grc_settings_file_path   : path to the GRC settings file.
      --chrom_key_file_path      : path to the chrom key file
      --kelch_reference_file_path: path to the kelch13 reference sequence file
      --codon_key_file_path      : path to the codon key file
      --drl_information_file_path: path to the drug resistance loci information file

    Settings:
      --results_dir     : output directory (default: output/)
      --panels_settings : path to panel_settings.csv
      --containers_dir  : path to a dir where the containers are located


    Additional options:
      --help    : prints this help message and exit
    
    Profiles:
      standard (default): run locally using singularity
      run_locally : run locally using what is available on the system environment (no containers)
   """.stripIndent()
}

def validateInputParams() {
  /*
  count errors on parameters which must be provided regardless of the workflow which will be executed
  
  returns
  -------
  <int> number of errors found
  */

  def error = 0
  def valid_execution_modes = ["cram", "fastq"]

  // check if execution mode is valid
  if (!valid_execution_modes.contains(params.execution_mode)){
    log.error("The execution mode provided (${params.execution_mode}) is not valid. valid modes = ${valid_execution_modes}")
    error += 1
  }

  // check if resources were provided
  error += __check_if_params_file_exist("grc_settings_file_path", params.grc_settings_file_path)
  error += __check_if_params_file_exist("panels_settings", params.panels_settings) 
  error += __check_if_params_file_exist("chrom_key_file_path", params.chrom_key_file_path) 
  error += __check_if_params_file_exist("codon_key_file_path", params.codon_key_file_path)
  error += __check_if_params_file_exist("drl_information_file_path", params.drl_information_file_path)

  if (params.no_kelch == false) {
    error += __check_if_params_file_exist("kelch_reference_file_path", params.kelch_reference_file_path)
  }
  
  // raise WARNING if debug params were set
  if (!params.DEBUG_takes_n_bams == null){
    log.warn("[DEBUG] takes_n_bams was set to ${params.DEBUG_takes_n_bams}")
  }

  if (!params.DEBUG_tile_limit == null){
    log.warn("[DEBUG] tile_limit was set to ${params.DEBUG_tile_limit}")
  }

  if (params.DEBUG_no_coi == true){
    log.warn("[DEBUG] no_coi was set to ${params.DEBUG_no_coi}")
  }
  // -------------------------------------------

  // check if output dir exists, if not create the default path
  if (params.results_dir){
    results_path = file(params.results_dir)
    if (!results_path.exists()){
      log.warn("${results_path} does not exists, the dir will be created")
      results_path.mkdir()
    }
  }

  if ((params.DEBUG_no_coi == false) && (params.mccoil_repopath != "/app/THEREALMcCOIL/")){
    mccoil_path = file(params.mccoil_repopath)
    if (mccoil_path.exists() == false){
      log.error("""
      The mccoil_repopath provided (${mccoil_path}) does not exists.
      This can happen if you do not use the containers provided or setup an invalid custom path.
      Please provide a valid custom installation path of the McCOIL library.
      """)
      error+=1
    }
  }

  if (error > 0) {
    log.error("Parameter errors were found, the pipeline will not run")
    exit 1
  }
} 

def __check_if_params_file_exist(param_name, param_value){
  // --- GRC SETTINGS ---
  def error = 0

  if (!(param_value==null)){
    param_file = file(param_value)
    if (!param_file.exists()){
      log.error("${param_file} does not exist")
      error +=1
    }
  }

  if (param_value==null){
    log.error("${param_name} must be provided")
    error +=1
  }
  // ----------------------
  return error
}

//
//-------------------- Validate Reference Panel ----------------------------------------
//
def addProjectDirAbsPathTo(inputString) {
    return inputString.replaceFirst("<ProjectDir>", "${projectDir}")
}

def replaceFileExtension(filePath, newExtension) {
    def lastIndex = filePath.lastIndexOf('.')
    if (lastIndex >= 0) {
        // Replace the old extension with the new one
        def newPath = filePath.substring(0, lastIndex) + '.' + newExtension
        return newPath
    } else {
        // If there is no existing extension, simply add the new extension
        return filePath + '.' + newExtension
    }
}


def validateReferenceFormat(row) {
    def errors = 0 

    // check if reference_file columns is a valid path
    if (row.reference_file.startsWith("<ProjectDir>")){
        reference_file = addProjectDirAbsPathTo(row.reference_file)
    } else {
        reference_file = row.reference_file
    }

    aligns_to_file = file(reference_file)

    if (!aligns_to_file.exists()){
        log.error("${reference_file} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if necessary index files for reference genome exist
    reference_index_file_list = [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]
    reference_index_file_list.each {extension  -> extension
        index_file = file("${reference_file}" + extension)
        if (!index_file.exists()) {
            log.error("${index_file} provided for ${row.panel_name} does not exist.")
            errors += 1
        }
    }

    // check if dictionary file for reference genome exists
    reference_dictionary_file = file(replaceFileExtension("${reference_file}", "dict"))
    if (!reference_dictionary_file.exists()){
        log.error("${reference_dictionary_file} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if snp_list columns is a valid path

    if (row.snp_list.startsWith("<ProjectDir>")){
        snp_list = addProjectDirAbsPathTo(row.snp_list)
    } else {
        snp_list = row.snp_list
    }

    snp_list_file = file(snp_list)
    if (!snp_list_file.exists()){
        log.error("${snp_list} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if design_file is a valid path
    if (row.design_file.startsWith("<ProjectDir>")){
        dsgn_path = addProjectDirAbsPathTo(row.design_file)
    } else {
        dsgn_path = row.design_file
    }

    maps_to_file = file(dsgn_path)
    if (!maps_to_file.exists()){
        log.error("${dsgn_path} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected at panel settings csv", errors))
        exit 1
    }
}

def validateReferencePanels(panels_settings) {
    // load panels settings content
    def panels_settings_ch = Channel.fromPath(panels_settings, checkIfExists: true)
                                | splitCsv(header: true, sep: ',') 

    // validate panels settings paths
    panels_settings_ch.map{row -> validateReferenceFormat(row)}
    // gen reference channel
    def reference_ch = panels_settings_ch 
                       | map { row ->
                         	// get absolute path for the pannel settings provided 
                            // on the repo. If not from the repo, load the path
                            // as set on the file
                            if (row.reference_file.startsWith("<ProjectDir>")){
                                ref_path = addProjectDirAbsPathTo(row.reference_file)
                            } else {
                            	ref_path = row.reference_file
                            }

                            if (row.snp_list.startsWith("<ProjectDir>")){
                                snp_path = addProjectDirAbsPathTo(row.snp_list)
                            } else {
                                snp_path = row.snp_list
                            }

                            tuple(row.panel_name, ref_path, snp_path)
                        }
    // gen annotations channel
    def annotations_ch = Channel.fromPath(panels_settings, checkIfExists: true)
                        | splitCsv(header: true, sep: ',')
                        | map { row ->
                            if (row.design_file.startsWith("<ProjectDir>")){
                                dsgn_path = file(addProjectDirAbsPathTo(row.design_file), checkIfExists: true)
                            } else {
                                dsgn_path = file(row.design_file, checkIfExists: true)
                            }
                            tuple( row.panel_name, dsgn_path )
                           }
    /* return [reference_ch, annotations_ch] */
    return reference_ch
}

//
//-----------Add reference--------------------
//

def groupPanelReference(input) {
	// Expect input of [primer_id, meta, input,fasta, snps]
    def (meta, input_files, fasta, snps) = input[1..4]
    reference_index_file_list = [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]
    genome = [
        fasta: fasta,
        fai: "${fasta}" + ".fai",
        amb: "${fasta}" + ".amb",
        ann: "${fasta}" + ".ann",
        bwt: "${fasta}" + ".bwt",
        pac: "${fasta}" + ".pac",
        sa:  "${fasta}" + ".sa",
        dict: "${fasta}" + ".dict",
    ]

    return [ meta + [ reference:genome, snps:snps ], input_files ]

}
