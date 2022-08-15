#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include { basecalls_conversion } from './modules/basecalls_conversion.nf'
include { decode_multiplexed_bam } from './modules/decode_multiplexed_bam.nf'
include { bam_find_adapter } from './modules/bam_find_adapter.nf'
include { bam_to_cram } from './modules/bam_to_cram.nf'
include { rename_cram_fls } from './modules/rename_cram_fls.nf'

process writeOutputManifest {

  publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

  input:
    tuple val(run_id), path(cram_files)
  
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

# get sample tags (nor proper ids yet, but will get the job done for now)
sample_tag_lst = [f.split("/")[-1].split(".")[0] for f in cram_files_lst]

# write manifest header

out_mnf.write("run_id,cram_fl,sample_tag\n")

# write manifest content
for idx in range(0, len(cram_files_lst)):
  cram_fl= cram_files_lst[idx]
  sample_tag = sample_tag_lst[idx]
  out_mnf.write(f"{run_id},{cram_dir}{cram_fl},{sample_tag}\n")
out_mnf.close()
/$
}

workflow bcl_to_cram {
    take:
        pre_process_input_ch
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
        bam_to_cram(bam_find_adapter.out.run_id,
                    bam_find_adapter.out.bam_adapter_file,
                    bam_find_adapter.out.bam_metrics_file)

        //cram_ch = bam_to_cram.out

        // rename samples to samplesheet provided names
        rename_cram_fls(bam_to_cram.out.run_id,
                        bam_to_cram.out.metrics_bam_file,
                        bam_to_cram.out.cram_fls
                        )
        cram_ch = rename_cram_fls.out

        // generate an output manifest
        writeOutputManifest(cram_ch)//cram_ch.run_id, cram_ch.cram_fls)
        manifest_out = writeOutputManifest.out

    emit:
        manifest_out
}
/*
// -------------------------- DOCUMENTATION -----------------------------------
[1] Mixed-sample MiSeq run sequencing data is first basecalled and then written to a BAM file.
[2] Resulting multiplexed BAM files are demultiplexed to separate reads based upon their sample of origin.
[3] BAM files are scanned for adapter sequence contamination, information about which is used to populate several auxiliary fields.
[4] BAM files are split by read group into CRAM files, which are emitted ready for further processing.

*/
