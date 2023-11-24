// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

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


def validatePanelSettings(row){
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
    if (!index_file.exists()){
        log.error("${index_file} provided for ${row.panel_name} does not exist.")
        errors += 1}
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

def parse_panel_settings(panels_settings) {
    def source_dir = ""
    // load panels settings content
    def panels_settings_ch = Channel.fromPath(panels_settings, checkIfExists: true)
                                | splitCsv(header: true, sep: ',') 

    // validate panels settings paths
    panels_settings_ch.map{row -> validatePanelSettings(row)}
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

                            tuple(ref_path,row.panel_name,snp_path)
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
                            tuple(
                                row.panel_name,
                                dsgn_path
                            )
                           }
    return [reference_ch, annotations_ch]
}
