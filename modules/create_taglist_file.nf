// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process create_taglist_file {

    input:
      val(study_name)
      path(manifest)

    output:
      path("${tag_list}")

    script:
      tag_list = "tag_file.tsv"

      """
      echo ${manifest}
      create_taglist_file.py -m $manifest --study $study_name
      """
}
