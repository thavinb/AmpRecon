#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process create_taglist_file {

    input:
      val(study_name)
      val(manifest)

    output:
      path("${tag_list}")

    script:
      tag_list = "tag_file.tsv"

      """
      echo ${manifest}
      create_taglist_file.py -m $manifest --study $study_name
      """
}
