/*
    | PARSE_PANNEL_SETTINGS |-----------------------------------------
    
    This workflow handles the pannel and resource bundle files used on 
    different steps of the pipeline. 

    Here, channels containing the pannels and its respective resource
    bundle files are obtained. 

    if a custom pannels settings is provided (as .csv) are 
    validated and parsed as channels to be emmited.

    if no custom pannels, channels will be created based on resource
    bundle available at the repo.

    ------------------------------------------------------------------
*/

// TODO: test if variation on "simple names" before the ."ext" breaks something on the pipeline
//       my guess is that there is some naming assumptions on some of the processes 

def validatePannelSettings(row){
    def errors = 0 
    
    //valid_pannel_names = ["PFA_GRC1_v1.0","PFA_GRC2_v1.0","PFA_Spec"]
    // check if pannel name is valid
    //if (!valid_pannel_names.contains(row.pannel_name)){
    //    log.error("${row.pannel_name} is not valid. Valid pannel names:\n ${valid_pannel_names}")
    //    errors += 1
    //}

    // check if align_to columns is a valid path
    align_to_path = file(row.aligns_to)
    if (!align_to_path.exists()){
        log.error("${row.aligns_to} provided for ${row.pannel_name} does not exist.")
        errors += 1
    }

    // check if maps_to_regions_of is a valid path
    maps_to_path = file(row.maps_to_regions_of)
    if (!maps_to_path.exists()){
        log.error("${row.maps_to_regions_of} provided for ${row.pannel_name} does not exist.")
        errors += 1
    }
    
    // TODO: check if all expected files are present on the resource bundle provided
    
    // count errors and kill nextflow if any had been found
    if (errors > 0) {
        log.error(String.format("%d errors detected at pannel settings csv", errors))
        exit 1
    }
}

workflow PARSE_PANNEL_SETTINGS {
    take:
        pannel_settings
        reference_dir // temporary
    main:
        // build reference channel from "aligns_to" column
        // if a csv was provided, mount channels for a given resource bundle
        if (!(pannels_settings==null)){
    
            Channel.fromPath(pannels_settings, checkIfExists: true)
                | splitCsv(header: true, sep: ',')
                | map { row ->
                        // validate row settings
                        validatePannelSettings(row)
                        // TODO: validate if expected files were provided
                        tuple(
                            file("${row.aligns_to}/*.fasta"),
                            row.pannel_name,
                            file("${row.aligns_to}/*.fasta.*")
                            )
                    }
                | set { reference_ch }
        }
        // if pannels settings csv is not provided, just use the files on the repo
        // TODO set a default behaviour based on a default pannel_settings.csv
        else {
            reference_ch = Channel.from(
                [file("${reference_dir}/grc1/*.fasta"), "PFA_GRC1_v1.0" , file("${reference_dir}/grc1/*.fasta.*")],
                [file("${reference_dir}/grc2/*.fasta"), "PFA_GRC2_v1.0", file("${reference_dir}/grc2/*.fasta.*")],
                [file("${reference_dir}/spec/*.fasta"), "PFA_Spec", file("${reference_dir}/spec/*.fasta.*")]
            )
        }

        // build pannel_anotations_files from "maps_to_regions_of"
        Channel.fromPath(pannels_settings, checkIfExists: true)
            | splitCsv(header: true, sep: ',')
            | map { row ->
                        // validate row settings
                        validatePannelSettings(row)
                        tuple(
                            row.pannel_name,
                            file("${row.maps_to_regions_of}")
                            )
                }
            | set { pannel_anotations_files }

    emit:
        reference_ch
        pannel_anotations_files
}
