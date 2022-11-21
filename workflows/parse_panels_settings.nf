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

    // check if align_to columns is a valid path
    aligns_to_path = "${source_dir}/${row.aligns_to}"
    aligns_to_file = file(aligns_to_path)
    if (!aligns_to_file.exists()){
        log.error("${aligns_to_path} provided for ${row.panel_name} does not exist.")
        errors += 1
    }

    // check if maps_to_regions_of is a valid path
    annotation_flpth = "${source_dir}/${row.maps_to_regions_of}"
    maps_to_file = file(annotation_flpth)
    if (!maps_to_file.exists()){
        log.error("${annotation_flpth} provided for ${row.panel_name} does not exist.")
        errors += 1
    }
    
    // TODO: check if all expected files are present on the resource bundle provided
    
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
                                file("${source_dir}/${row.aligns_to}/*.fasta"),
                                row.panel_name,
                                file("${source_dir}/${row.aligns_to}/*.fasta.*")
                            )
                            }

        // build panel_anotations_files from "maps_to_regions_of"
        annotations_ch = Channel.fromPath(panels_settings, checkIfExists: true)
                          | splitCsv(header: true, sep: ',')
                          | map { row ->
                                // validate row settings
                                validatePanelSettings(row, source_dir)
                                tuple(
                                    row.panel_name,
                                    file("${source_dir}/${row.maps_to_regions_of}")
                                )
                           }

    emit:
        reference_ch
        annotations_ch
}
