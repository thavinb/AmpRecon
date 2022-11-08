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

    // check if reference_index_files is a valid path
    reference_index = "${source_dir}/${row.reference_index_file_basename}{.fai,.amb,.ann,.bwt,.pac,.sa}"
    reference_index_file_list = file(reference_index)
    reference_index_file_list.each { file -> if (!file.exists()){
        log.error("${file} provided for ${row.panel_name} does not exist.")
        errors += 1}
    }

    // check if reference_dictionary_file is a valid path
    reference_dictionary = "${source_dir}/${row.reference_dictionary_file}"
    reference_dictionary_file = file(reference_dictionary)
    if (!reference_dictionary_file.exists()){
        log.error("${reference_dictionary} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if annotation_vcf_file columns is a valid path
    annotation_vcf = "${source_dir}/${row.annotation_vcf_file}"
    annotation_vcf_file = file(annotation_vcf)
    if (!annotation_vcf_file.exists()){
        log.error("${annotation_vcf} provided for ${row.panel_name} does not exist.")
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
                                file("${source_dir}/${row.reference_file}"),
                                row.panel_name,
                                file("${source_dir}/${row.reference_index_file_basename}{.fai,.amb,.ann,.bwt,.pac,.sa}"),
                                file("${source_dir}/${row.reference_dictionary_file}"),
                                file("${source_dir}/${row.ploidy_file}"),
                                file("${source_dir}/${row.annotation_vcf_file}"),
				file("${source_dir}/${row.snp_list}")
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
//reference_ch.first().view()
    emit:
        reference_ch
        annotations_ch // tuple ([fasta], panel_name, [fasta_idx_files], [dictionary_file], [ploidy_file], [annotation_vcf_file], [snp_list])
}
