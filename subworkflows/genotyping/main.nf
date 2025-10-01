// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include { BCFTOOLS_MPILEUP } from '../..//modules/bcftools_mpileup/main'
include { BCFTOOLS_CALL    } from '../../modules/bcftools_call/main'
include { BCFTOOLS_FILTER  } from '../../modules/bcftools_filter/main'
include { BCFTOOLS_STATS   } from '../../modules/bcftools_stats/main'

workflow GENOTYPING {

  take:
    bam_ch // tuple(meta, bam, bai)

  main:
    
    ch_versions = Channel.empty()
	ch_multiqc_files = Channel.empty()

    // MPILEUP
    BCFTOOLS_MPILEUP(bam_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    // call SNP sites 
    BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.bcf)
    ch_versions = ch_versions.mix(BCFTOOLS_CALL.out.versions.first())

    // filter genetic variants
    BCFTOOLS_FILTER(BCFTOOLS_CALL.out.bcf)
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

	// Gettiing variant call statistics
	BCFTOOLS_STATS(BCFTOOLS_FILTER.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())
	ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{it[1]})

  emit:
    vcf = BCFTOOLS_FILTER.out.vcf // tuple(meta, vcf, tbi)
	mqc = ch_multiqc_files
    versions = ch_versions

}

