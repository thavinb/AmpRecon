#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

        // container should be included later

        download_bcl_from_s3("21045") //gets test bcl directory from s3
        write_manifest("21045")
}

process download_bcl_from_s3 {
        publishDir "${params.test_data_dir}", mode: 'move'

	//download and unzip bcl related files for testing purposes

	input:
            val(bcl_id)

	output:
            path("${bcl_id}")

	script:
            """
            s3cmd get s3://amplicon-test-data/${bcl_id}.tar.gz
            tar -xvzf ${bcl_id}.tar.gz
            """
}

process write_manifest {
        publishDir "${params.test_data_dir}", mode: 'copy'

        input:
            val(run_id)

        output:
            path("manifest.csv")

// The $/ ... /$ is necessary to avoid nextflow to read "\n" correctly
$/
#!/usr/bin/python3

# setup inputs
run_id = "${run_id}"
out_mnf = open("manifest.csv", "w")
bcl_dir=f"${params.test_data_dir}/"
lane = 1
study_name = "test"
read_group = "${run_id}_1"
library = "lib"

# write manifest header

out_mnf.write("run_id,bcl_dir_path,lane,study_name,read_group,library\n")

# write manifest content
out_mnf.write(f"{run_id},{bcl_dir}{run_id},{lane},{study_name},{read_group},{library}\n")
out_mnf.close()
/$

}
