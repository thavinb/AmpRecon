#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { basecalls_conversion } from '../modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from '../modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from '../modules/bam_find_adapter.nf'
include { bam_to_cram } from '../modules/bam_to_cram.nf'

process writeOutputManifest {
  publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

input:
   tuple val(run_id), path(cram_files)
   // TODO create a python box container

output:
   tuple val(run_id), path("${run_id}_out1.1_mnf.csv")
// The $/ ... /$ is necessary to avoid nextflow to read "\n" correctly
$/
#!/usr/bin/python3

# setup inputs
run_id = "${run_id}"
cram_files_lst = "${cram_files}".split(" ")
out_mnf = open("${run_id}_out1.1_mnf.csv", "w")
cram_dir=f"${params.results_dir}{run_id}/"

# write manifest header

out_mnf.write("run_id,cram_fl\n")

# write manifest content
for cram_fl in cram_files_lst:
  out_mnf.write(f"{run_id},{cram_dir}{cram_fl}\n")
out_mnf.close()
/$
}

workflow bcl_to_cram {
    take:
        // pre_process_input_ch
        pre_process_input_ch
        //bcl_dir
        //lane
        //read_group
        //study_name
        //barcode
    main:
        // convert basecalls
        basecalls_conversion(pre_process_input_ch)
        decode_In_ch = pre_process_input_ch.join(basecalls_conversion.out)
        // decode multiplexed bam file
        decode_multiplexed_bam(decode_In_ch)
        find_adapter_In_ch = decode_multiplexed_bam.out

        // find adapter contamination in bam
        bam_find_adapter(find_adapter_In_ch)
        // split bam by read group into cram
        bam_to_cram(bam_find_adapter.out)
        cram_ch = bam_to_cram.out
        writeOutputManifest(cram_ch)
        manifest_out = writeOutputManifest.out
        // generate an output manifest

    emit:
        manifest_out
}
