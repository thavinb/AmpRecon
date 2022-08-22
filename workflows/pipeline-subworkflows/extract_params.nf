include { MAKE_SAMPLESHEET_MANIFEST } from '../../modules/make_samplesheet_manifest.nf'
include { VALIDATE_SAMPLESHEET_MANIFEST } from '../../modules/samplesheet_manifest_validation.nf'
include { GET_TAG_LIST_FILE } from '../../modules/manifest2tag.nf'

workflow EXTRACT_PARAMS {

    main:
    Channel.fromPath(params.input_params_csv)
           .splitCsv(header: true)
           .map { it -> tuple ( it.run_id,
                                it.bcl_dir_path,
                                it.lane,
                                it.study_name,
                                it.read_group )
          }
          .set { output_ch }

    barcode_gen = output_ch.map { it -> tuple( it[0], it[1] ) }
    study_name = output_ch.map { it -> it[3] }

    //get barcodes
    barcodes_file = MAKE_SAMPLESHEET_MANIFEST(barcode_gen).manifest_file

    //validate manifest and get tags
    VALIDATE_SAMPLESHEET_MANIFEST(barcodes_file)
    GET_TAG_LIST_FILE(VALIDATE_SAMPLESHEET_MANIFEST.out, "library", "sample", study_name)
    tag_list_file = GET_TAG_LIST_FILE.out.taglist_file

    Channel.fromPath(params.input_params_csv)
           .splitCsv(header : true)
           .map { it -> tuple( it.run_id, it.reference_fasta ) }
           .set { reference_information }


    emit:
    output_ch
    barcodes_file
    tag_list_file
    reference_information


}
