workflow ASPIS_0_1_0_Stage1_Step2 {

  File Manifest
  Array[File] InputCrams
  String RefPath

  scatter (cram in InputCrams) {

    call get_reference as ReferenceFile {
      input:
        Manifest = Manifest,
        InputCram = cram,
        ReferencePath = RefPath
    }

    call task_1_collate {
      input:
        InputCram = cram
    }

    call task_2_bamrest {
      input:
        InputCram = task_1_collate.OutputCram
    }

    call task_3_adapter_clip {
      input:
        InputCram = task_2_bamrest.OutputCram
    }

    call task_4_cram_2_fastq {
      input:
        InputCram = task_3_adapter_clip.OutputCram
    }

    call task_5_bwa_align {
      input:
        InputFastq = task_4_cram_2_fastq.OutputFastq,
        Reference = ReferenceFile.Path
    }

    call task_6_scramble {
      input:
        InputSam = task_5_bwa_align.OutputSam
    }

    call build_header as Task7Results {
      input:
        Task3InputCram = task_3_adapter_clip.OutputCram,
        Task6InputBam = task_6_scramble.OutputBam,
        Reference = ReferenceFile.Path
    }

    call task_7_reheader {
      input:
        InputBam = task_6_scramble.OutputBam,
        Header = Task7Results.Header
    }

    call task_8_bam_split {
      input:
        InputBam = task_7_reheader.OutputBam
    }

    call task_9_bam12uxmerger {
      input:
        InputBam = task_8_bam_split.OutputBam,
        Task3InputBam =  task_3_adapter_clip.OutputCram
    }

    call task_10_alignment_filter {
      input:
        InputBam = task_9_bam12uxmerger.OutputBam,
        Task1InputCram = task_1_collate.OutputCram
    }

    call task_11_bamsort {
      input:
        InputBam = task_10_alignment_filter.OutputBam
    }

    call task_12_scramble2 {
      input:
        InputBam = task_11_bamsort.OutputBam,
        Reference = ReferenceFile.Path
    }
  }

  output {
    #A list of crams
    Array[File] BamList = task_12_scramble2.OutputBam
    Array[String] ReferenceList = ReferenceFile.Path
  }

}

# --- Helper Tasks ---------------
task get_reference {
  File Manifest
  File InputCram
  String ReferencePath

  command <<<
  python3 <<CODE
  #load manifest
  from aspismanifest import ASPISManifest, ASPISManifestParser
  import json, os
  parser = ASPISManifestParser()
  man = parser.parse("${Manifest}")
  _refPath = "${ReferencePath}"

  if man[0] is None:
   os.exit(-1)
  man = man[0]

  #convert input cram into
  cramPath = "${InputCram}"
  index = cramPath[cramPath.rfind("#") + 1: -5]

  refPath = ""
  for sample in man.data:
   if index == sample["index"]:
    if sample["ref"] == "PFA_GRC1_v1.0":
     refPath = os.path.join(_refPath,"grc1","Pf_GRC1v1.0.fasta")
    elif sample["ref"] == "PFA_Spec":
     refPath = os.path.join(_refPath,"spec","Spec_v1.0.fasta")
    elif sample["ref"] == "PFA_GRC2_v1.0":
     refPath = os.path.join(_refPath,"grc2","Pf_GRC2v1.0.fasta")
    break
  print(refPath)

  CODE
  >>>

  output {
    String Path = read_string(stdout())
  }

  runtime {
    docker:"aspis-tools:latest"
  }
}

# --------------------------------
task build_header {
  File Task3InputCram
  File Task6InputBam
  String Reference

  command <<<
  python3 <<CODE
  import sys
  import time
  import os
  import shutil
  import subprocess
  import re

  def querySamtools(target):
   result = subprocess.run(["samtools","view","-H",target],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   if result.returncode == 0:
    return result.stdout.decode("UTF-8")
   else:
    return None
  #-------------------------------------------------

  def getHeaderValue(line,target):
   if target in line:
    start = line.index(target)
    end   = line.index("\t",start)
    return line[start:end]
   return ""
  #-------------------------------------------------
  regex = re.compile("@SQ.*\tSN:([^\t]+)")

  old = querySamtools("${Task3InputCram}")
  new = querySamtools("${Task6InputBam}")

  lastPGID = ""
  filteredOldArray = []
  for l in old.split("\n"):
   filteredOldArray.append(l)
   if l.startswith("@PG") and "ID:" in l:
    lastPGID = getHeaderValue(l,"ID:").replace("ID","\tPP")

  old = "\n".join(filteredOldArray)

  lookup = {}
  with open("${Reference}.dict","r") as fd:
   for line in fd:
    if "\t" in line:
     refBits = line.split("\t")
     lookup[refBits[1]] = line.strip()

   filteredNewArray = []
   newPGs = []
   newlines = new.split("\n")
   for line in newlines:
    if regex.search(line):
     refBits = line.split("\t")
     if refBits[1] in lookup:
      filteredNewArray.append(lookup[refBits[1]])
    else:
     if "@PG" in line:
      if "PP:" not in line:
       newPGs.append(line + lastPGID)
      else:
       newPGs.append(line)
     else:
      filteredNewArray.append(line)

  new = "\n".join(filteredNewArray)

  sys.stdout.write(old.strip())
  sys.stdout.write("\n")
  sys.stdout.write("\n".join(reversed(newPGs)).strip())
  sys.stdout.write("\n")
  sys.stdout.write(new.strip())
  sys.stdout.write("\n")

  CODE
  >>>

  output {
    String Header = read_string(stdout())
  }

  runtime {
    docker:"samtools:latest"
    docker_args:"-v /home/ubuntu/inputs/references:/home/ubuntu/inputs/references"
  }
}

# --- Tasks ----------------------
task task_1_collate {
  File InputCram #$1/step1-4/$2_1#$3.cram

  command <<<
  bamcollate2 collate=1 level=0 \
  inputformat=cram \
  filename=${InputCram} \
  O=${basename(InputCram)}
  >>>

  output {
    File OutputCram = basename(InputCram)
  }

  runtime {
    docker:"biobambam2:latest"
  }
}

# --------------------------------
task task_2_bamrest {
  File InputCram #$1/step2-1/$2_1#$3.cram

  command <<<
  bamreset \
  resetaux=0 \
  level=0 \
  verbose=0 \
  inputformat=cram \
  < ${InputCram} \
  > ${basename(InputCram)}
  >>>

  output {
    File OutputCram = basename(InputCram)
  }

  runtime {
    docker:"biobambam2:latest"
  }
}

# --------------------------------
task task_3_adapter_clip {
  File InputCram #$1/step2-2/$2_1#$3.cram

  command <<<
  bamadapterclip verbose=0, level=0, \
  <  ${InputCram} \
  > ${basename(InputCram)}
  >>>

  output {
    File OutputCram = basename(InputCram)
  }

  runtime {
    docker:"biobambam2:latest"
  }
}

# --------------------------------
task task_4_cram_2_fastq {
  File InputCram #$1/step2-3/$2_1#$3.cram
  String outputName = sub(basename(InputCram),"cram$","fastq") #$1/step2-4/$2_1#$3.fastq

  command <<<
  bamtofastq \
  < ${InputCram} \
  > ${outputName}
  >>>

  output {
    File OutputFastq = outputName
  }

  runtime {
    docker:"biobambam2:latest"
  }
}

# --------------------------------
task task_5_bwa_align {
  File InputFastq #$1/step2-4/$2_1#$3.fastq
  String Reference #$4
  String outputName =  sub(basename(InputFastq),"fastq$","sam") #$1/step2-5/$2_1#$3.sam

  command <<<
  bwa mem \
  -t 24 \
  -p \
  -Y \
  -K 100000000 \
  ${Reference} \
  ${InputFastq} \
  > ${outputName}
  >>>

  output {
    File OutputSam = outputName
  }

  runtime {
    docker:"bwa:latest"
    docker_args:"-v /home/ubuntu/inputs/references:/home/ubuntu/inputs/references"
  }
}

# --------------------------------
task task_6_scramble {
  File InputSam #$1/step2-5/$2_1#$3.sam
  String outputName = sub(basename(InputSam),"sam$","bam") #$1/step2-6/$2_1#$3.bam

  command <<<
  scramble \
  -0 \
  -I sam \
  -O bam \
  < ${InputSam} \
  > ${outputName}
  >>>

  runtime {
    docker:"staden:latest"
  }

  output {
    File OutputBam = outputName
  }
}

# --------------------------------
task task_7_reheader {
  File InputBam
  String Header
  String outputName = basename(InputBam)

  command <<<
  echo "${Header}" | samtools reheader - \
  ${InputBam} \
  > ${outputName}
  >>>

  output {
    File OutputBam = outputName
  }

  runtime {
    docker:"samtools:latest"
  }
}

# --------------------------------
task task_8_bam_split {
  File InputBam #$1/step2-7/$2_1#$3.bam
  String outputName = basename(InputBam) #$1/step2-8/$2_1#$3.bam

  command <<<
  bam12split verbose=0 level=0 \
  < ${InputBam} \
  > ${outputName}
  >>>

  output {
    File OutputBam = outputName
  }

  runtime {
    docker:"biobambam2:latest"
    docker_args:"-v /home/ubuntu/inputs/references:/home/ubuntu/inputs/references"
  }
}

# --------------------------------
task task_9_bam12uxmerger {
  File InputBam #$1/step2-8/$2_1#$3.bam
  File Task3InputBam #$1/step2-3/$2_1#$3.cram
  String outputName = basename(InputBam) #$1/step2-9/$2_1#$3.bam

  command <<<
  bam12auxmerge \
  level=0 \
  rankstrip=1 \
  ranksplit=0 \
  zztoname=0 \
  clipreinsert=1 \
  ${Task3InputBam} \
  < ${InputBam} \
  > ${outputName}
  >>>

  output {
    File OutputBam = outputName
  }

  runtime {
    docker:"biobambam2:latest"
  }
}

# --------------------------------
task task_10_alignment_filter {
  File InputBam #$1/step2-9/$2_1#$3.bam
  File Task1InputCram #$1/step2-1/$2_1#$3.cram
  String outputName = basename(InputBam) #$1/step2-10/$2_1#$3.bam
  String splitName = sub(outputName,"\\.bam$","_split.bam") #$1/step2-10/$2_1#$3-phix.bam
  String metricsName = sub(outputName,"\\.bam$","_bam_alignment_filter_metrics.json") #$1/step2-10/$2_1#$3_bam_alignment_filter_metrics.json

  command <<<
  bambi select \
  --compression-level=0 \
  --input ${Task1InputCram},${InputBam} \
  --output ${splitName},${outputName} \
  -m ${metricsName}
  >>>

  output {
    File OutputBam = outputName
    File BamAlignmentFilterMetrics = metricsName
  }

  runtime {
    docker:"bambi:latest"
  }
}

# --------------------------------
task task_11_bamsort {
  File InputBam #$1/step2-10/$2_1#$3.bam
  String outputName = basename(InputBam) #$1/step2-11/$2_1#$3.bam

  command <<<
  bamsort threads=24 SO=coordinate level=0 fixmate=1 \
  addupmarksupport=1 tmpfile=9_bamsort.tmp \
  < ${InputBam} \
  > ${outputName}
  >>>

  output {
    File OutputBam = outputName
  }

  runtime {
    docker:"biobambam2:latest"
  }
}

# --------------------------------
task task_12_scramble2 {
  String Reference #$4
  File InputBam #$1/step2-11/$2_1#$3.bam
  String outputName = basename(InputBam) #$1/step2-12/$2_1#$3.bam

  command <<<
  scramble -t 11 -7 -I bam -O cram \
  -r ${Reference} \
  -e \
  < ${InputBam} \
  > ${outputName}
  >>>

  output {
    File OutputBam = outputName
  }

  runtime {
    docker:"staden:latest"
    docker_args:"-v /home/ubuntu/inputs/references:/home/ubuntu/inputs/references"
  }
}
