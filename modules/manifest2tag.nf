#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process get_taglist_file {

    input:
      val(run_id)
      val(study_name)
      val(library)
      val(manifest)

    output:
      path("${tag_list}")

    script:
      tag_list = "tag_file.tsv"

      """
      echo ${manifest}
      manifest2tag.py -m $manifest -l $library --sample $run_id --study $study_name
      """
}
