include { CORE_PIPELINE_REPLICA } from './pipeline-subworkflows/core_pipeline_replica.nf'
include { EXTRACT_PARAMS } from './pipeline-subworkflows/extract_params.nf'


workflow IN_COUNTRY {

    main:
    EXTRACT_PARAMS_INCOUNTRY()
    CORE_PIPELINE_REPLICA(EXTRACT_PARAMS.out)

}
