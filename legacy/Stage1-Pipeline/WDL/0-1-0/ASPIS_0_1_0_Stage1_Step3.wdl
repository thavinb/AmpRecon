workflow ASPIS_0_1_0_Stage1_Step3 {

    Array[File] InputBams
    Array[String] References

    scatter (index in range(length(InputBams)))
    {
      call task_1_mpileup {
        input:
          InputBam = InputBams[index],
          Reference = References[index]
      }

      call task_2_call {
        input:
          InputBCF = task_1_mpileup.OutputBCF,
          Reference = References[index]
      }

      call task_3_filter_depth {
        input:
          InputBCF = task_2_call.OutputBCF
      }

      call task_4_filter_quality {
        input:
          InputVCF = task_3_filter_depth.OutputVCF
      }

      output {
        Array[File] VCFList = task_4_filter_quality.OutputVCF
      }
    }

}

# --------------------------------
task task_1_mpileup {

  File InputBam
  String outputName = sub(basename(InputBam),"bam$","bcf")
  String Reference
  String annotation = sub(Reference,"fasta$","annotation.vcf")

  command <<<
  bcftools mpileup \
  --min-BQ 20 \
  --annotate FORMAT/AD,FORMAT/DP \
  --max-depth 50000 \
  --targets-file ${annotation} \
  --fasta-ref ${Reference} \
  --output-type u \
  ${InputBam} \
  > ${outputName}
  >>>

  runtime {
    docker:"bcftools:latest"
    docker_args:"-v /home/ubuntu/inputs/references:/home/ubuntu/inputs/references"
  }

  output {
    File OutputBCF = outputName
  }
}

# --------------------------------
task task_2_call {

  File InputBCF
  String outputName = basename(InputBCF)
  String Reference
  String ploidy = sub(Reference,"fasta$","ploidy")

  command <<<
  bcftools call \
  --multiallelic-caller \
  --keep-alts \
  --skip-variants indels \
  --ploidy-file ${ploidy} \
  --output-type u \
  < ${InputBCF} \
  > ${outputName}
  >>>

  runtime {
    docker:"bcftools:latest"
    docker_args:"-v /home/ubuntu/inputs/references:/home/ubuntu/inputs/references"
  }

  output {
    File OutputBCF = outputName
  }
}

# --------------------------------
task task_3_filter_depth {

  File InputBCF
  String outputName = sub(basename(InputBCF),"bcf$","vcf")

  command <<<
  bcftools filter \
  --mode + \
  --soft-filter LowDepth \
  --exclude FORMAT/DP\<8 \
  --output-type v \
  < ${InputBCF} \
  > ${outputName}
  >>>

  runtime {
    docker:"bcftools:latest"
  }

  output {
    File OutputVCF = outputName
  }
}

# --------------------------------
task task_4_filter_quality {

  File InputVCF
  String outputName = basename(InputVCF)

  command <<<
  bcftools filter \
  --mode + \
  --soft-filter LowQual \
  --exclude '%QUAL<15 || MQ<20' \
  --output-type v \
  < ${InputVCF} \
  > ${outputName}
  >>>

  runtime {
    docker:"bcftools:latest"
  }

  output {
    File OutputVCF = outputName
  }
}
