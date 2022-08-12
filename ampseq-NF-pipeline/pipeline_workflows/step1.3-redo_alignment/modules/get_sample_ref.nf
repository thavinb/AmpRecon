process get_sample_ref {

    //publishDir "${params.results_dir}/${run_id}", mode: 'copy', overwrite: true

    input:
        val(run_id)
        val(sample_tag)
        //path(manifest) // e.g. bam
        path(cram_fl)
    output:
        tuple val(run_id), val(cram_fl), val(sample_tag), path("*.fasta"), path("*.fasta.*"), path("*.fasta.fai"), path("*.fasta.dict")
        //val(run_id), emit: run_id
        //val(cram_fl), emit: cram_fl
        //val(sample_tag), emit: sample_tag
        //path("*.fasta"), emit: ref_fasta_file
        // this will include fai, needs to improve this in the future
        //path("*.fasta.*"), emit: ref_bwa_index_fls
        //path("*.fasta.fai"), emit: ref_fasta_idx_fl
        //path("*.fasta.dict"), emit: ref_fasta_dct
script:
// GAMBIARRA ALERT ------------------------------------------------------------
manifest = "${params.results_dir}/${run_id}/${run_id}_manifest.csv"
// ----------------------------------------------------------------------------
ref_database = "${projectDir}/references/"
"""
python3 ${projectDir}/pipeline_workflows/step1.2a-cram-to-bam/modules/scripts/getSamplesRef.py \
              -m ${manifest} -cfl ${cram_fl} -rdb ${ref_database}
"""
}
