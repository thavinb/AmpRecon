## readManifest 0-1-0 - Reads a ASPIS manifest and generates a manifest object

workflow loadManifest {
  File Manifest

  call parseAndSplit {
    input:
      Manifest = Manifest
  }

  call createTaglist {
    input:
      Manifest = Manifest,
      Batch = parseAndSplit.Header["batch"]
  }

  output {
    File Taglist = createTaglist.Taglist
    Array[String] Samples = parseAndSplit.Samples
    Map[String,String] Header  = parseAndSplit.Header
  }
}

#Use manaifest parser to load manifest and then
task parseAndSplit {
  File Manifest

  command <<<
  python3 <<CODE
  from aspismanifest import ASPISManifest, ASPISManifestParser
  import json, os
  parser = ASPISManifestParser()
  man = parser.parse("${Manifest}")

  if man[0] is None:
   os.exit(-1)

  man = man[0]

  headers = {}
  headers["lane"] = man.lane
  headers["pipeline"]= "{}_{}".format(man.pipeline,man.pipelineVersion)
  headers["batch"]= man.batch
  headers["study"]= "ASPIS_0_1_0"

  with open('headers.json', 'w', encoding='utf-8') as f:
   json.dump(headers, f, ensure_ascii=False, indent=4)

  with open('samples.txt', 'w', encoding='utf-8') as f:
   control = 1
   for r in man.data:
    sim = r["sims"]
    if sim.lower() == "control":
     sim = "Control-{}".format(control)
     control += 1
    f.write("{}\n".format(sim))

  CODE
  >>>

  output {
    Array[String] Samples = read_lines('samples.txt')
    Map[String,String] Header  = read_json('headers.json')
  }

  runtime {
    docker:"aspis-tools:latest"
  }
}

task createTaglist {
  File Manifest
  String Batch

  command <<<
  aspis-tools manifest2taglist ${Manifest}
  >>>

  output {
    File Taglist = "${Batch}.taglist"
  }

  runtime {
    docker:"aspis-tools:latest"
  }
}
