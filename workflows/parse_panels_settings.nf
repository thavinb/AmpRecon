/*
    | PARSE_PANEL_SETTINGS |-----------------------------------------
    
    This workflow handles the panel and resource bundle files used on 
    different steps of the pipeline. 

    Here, channels containing the panels and its respective resource
    bundle files are obtained. 

    if a custom panels settings is provided (as .csv) are 
    validated and parsed as channels to be emmited.

    if no custom panels, channels will be created based on resource
    bundle available at the repo.

    ------------------------------------------------------------------
*/

// TODO: test if variation on "simple names" before the ."ext" breaks something on the pipeline
//       my guess is that there is some naming assumptions on some of the processes 

def validatePanelSettings(row, source_dir){
    def errors = 0 

    // check if reference_file columns is a valid path
    reference_file = "${source_dir}/${row.reference_file}"
    aligns_to_file = file(reference_file)
    if (!aligns_to_file.exists()){
        log.error("${reference_file} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if necessary index files for reference genome exist
    reference_index_file_list = [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]
    reference_index_file_list.each {extension  -> extension
    index_file = file("${source_dir}/${row.reference_file}" + extension)
    if (!index_file.exists()){
        log.error("${index_file} provided for ${row.panel_name} does not exist.")
        errors += 1}
    }

    // check if dictionary file for reference genome exists
    reference_dictionary = "${source_dir}/" + file("${row.reference_file}").parent + "/" + file("${row.reference_file}").baseName+".dict"
    reference_dictionary_file = file(reference_dictionary)
    if (!reference_dictionary_file.exists()){
        log.error("${reference_dictionary} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if snp_list columns is a valid path
    snp_list = "${source_dir}/${row.snp_list}"
    snp_list_file = file(snp_list)
    if (!snp_list_file.exists()){
        log.error("${snp_list} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if design_file is a valid path
    annotation_flpth = "${source_dir}/${row.design_file}"
    maps_to_file = file(annotation_flpth)
    if (!maps_to_file.exists()){
        log.error("${annotation_flpth} provided for ${row.panel_name} does not exist.")
        errors += 1
    }
    
    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected at panel settings csv", errors))
        exit 1
    }
}

workflow PARSE_PANEL_SETTINGS {
    take:
        panels_settings
    
    main:
        // if a csv was provided, no need to add source_dir to rows
        if (!(panels_settings==null)){
            source_dir = ""
        }
        // if no panels_settings csv is provided, use the one at the repo
        else {
            source_dir = "${projectDir}" // required to get the right path of resources at repo
            panels_settings = "${source_dir}/panels_resources/panels_settings.csv"
        }
        
        // build reference channel from "aligns_to"
        reference_ch = Channel.fromPath(panels_settings, checkIfExists: true)
                        | splitCsv(header: true, sep: ',')
                        | map { row ->
                            // validate row settings
                    
                            validatePanelSettings(row, source_dir)
                            // TODO: validate if expected files were provided
                            tuple(
                                "${source_dir}/${row.reference_file}",
                                row.panel_name,
				"${source_dir}/${row.snp_list}"
                            )
                            }

        // build panel_anotations_files from "design_file"
        annotations_ch = Channel.fromPath(panels_settings, checkIfExists: true)
                          | splitCsv(header: true, sep: ',')
                          | map { row ->
                                // validate row settings
                                validatePanelSettings(row, source_dir)
                                tuple(
                                    row.panel_name,
                                    file("${source_dir}/$row.design_file")
                                )
                           }
    emit:
        reference_ch // tuple(reference_file, panel_name, snp_list)
        annotations_ch // tuple(panel_name, design_file)
}
