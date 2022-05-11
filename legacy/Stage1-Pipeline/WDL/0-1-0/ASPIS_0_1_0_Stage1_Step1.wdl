## Stage 1 BCL To VCF Step 1 - Convert BCL into a set of BAMS

import "readManifest.wdl" as readManifest

# Inputs - Mainfest, bcl path
# Manifest becomes - taglist
workflow ASPIS_0_1_0_Stage1_Step1 {
  File BCLFolder
  File Manifest

  call readManifest.loadManifest as ManifestData {
    input:
      Manifest = Manifest
  }

  call BCLtoBAM {
    input:
      BCLFolder = BCLFolder,
      Lane = ManifestData.Header["lane"],
      PlatformUnit = ManifestData.Header["pipeline"],
      Study = ManifestData.Header["study"],
      Samples = ManifestData.Samples,
      Batch = ManifestData.Header["batch"]
  }

  call BambiDecode {
    input:
      Taglist = ManifestData.Taglist,
      InputBam = BCLtoBAM.Result
  }

  call BamAdapterFind {
    input:
      InputBam = BambiDecode.Result,
  }

  call SamToolsSplit {
    input:
      InputBam = BamAdapterFind.Result
  }

  output {
    Array[File] CramList = SamToolsSplit.CramList
  }

}

#Task 1 Convert BCL to BAM -----------------------------------------------------
task BCLtoBAM {
  String BCLFolder
  Int Lane
  String PlatformUnit
  String Study
  Array[String] Samples
  String Batch

  command <<<
    bambi i2b --intensity-dir=/blc/Data/Intensities \
     --basecalls-dir=/blc/Data/Intensities/BaseCalls \
     --lane=${Lane} \
     --platform-unit=${PlatformUnit} \
     --read-group-id=${Batch}_${Lane} \
     --study-name=${Study} \
     --sample-alias=${sep="," Samples} \
     --threads=8 \
     --output-file=step1-1.bam \
     --compression-level=0
  >>>

 runtime {
   docker:"bambi:latest"
   docker_args:"-v ${BCLFolder}:/blc"

 }

 output {
   File Result = "step1-1.bam"
 }

}

#Task 2 Bambi decode -----------------------------------------------------------
task BambiDecode {
  File Taglist
  File InputBam

  command <<<
  bambi decode --output=step1-2.bam \
 --compression-level=0 \
 --metrics-file=matrix.txt \
 --barcode-file=${Taglist} \
 ${InputBam}
 >>>

 runtime {
   docker:"bambi:latest"
 }

 output {
   File Result = "step1-2.bam"
 }
}

#Task 3 find adapters ----------------------------------------------------------
task BamAdapterFind {

  File InputBam

  command <<<
  bamadapterfind level=0 < ${InputBam}  > step1-3.bam
  >>>

  runtime {
    docker:"biobambam2:latest"
  }

  output {
    File Result = "step1-3.bam"
  }
}

#Task 4 samtools - split -------------------------------------------------------
task SamToolsSplit {

  File InputBam

  command <<<
  samtools split --threads 4 -v \
  --output-fmt cram,no_ref=1 \
  -f %!.cram \
  ${InputBam}

  ls -1 | grep .cram >> "cramlist.txt"
  >>>

  runtime {
    docker:"samtools:latest"
  }

  output {
    Array[File] CramList = read_lines("cramlist.txt")
  }

}
