import os
import sys
import json
import shutil
import subprocess

from aspismanifest import ASPISManifest, ASPISManifestParser
from aspismanifesttools import manifest2taglist

#------------------------------------------------- 
def executeStep(_target,_args,_execPath):
 print("Target: {} ExecPath: {}".format(_target,_execPath))
 args = ["bash",_target +".sh"]
 args.extend(_args)
 # test ---------------------------------------------------
 container_path ="/home/ubuntu/Projects/containers/Bambi/bambi.sif"
 args = ["bash", "singularity", "exec", container_path, "/bin/sh", _target+".sh"]
 args.extend(_args)
 print(" ".join(args))
 print(stop)
 # --------------------------------------------------------
 result = subprocess.run(" ".join(args),stderr=subprocess.PIPE,cwd=_execPath,shell=True)
 if result.returncode != 0:
  exitWithError("Process failed with error {}\n{}".format(result.returncode,result.stderr.decode("UTF-8")))
  return False
 return True

#------------------------------------------------- 

def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()

#------------------------------------------------- 

def runStep1(_manifest,_executionPath,_bclPath, _savePath):
 # create save paths
 if not os.path.isdir(_savePath):
  return (False,"Save path does not exist")
 
 call4SavePath = os.path.join(_savePath,"step1-4")
 if not os.path.exists(call4SavePath):
  os.makedirs(call4SavePath)
  
 # create taglist
 taglist = manifest2taglist(_manifest,_savePath,_manifest.batch)
 print("Taglist saved to {}".format(taglist))
 
 sampleNames = []
 controlIndex = 1
 for d in _manifest.data:
  name = d["sims"]
  if name.lower() == "control":
   name = "Control-{}".format(controlIndex)
  controlIndex += 1
  sampleNames.append(name)
 
 # call 1 - bcl root path, platform unit, lane, batch, study, comma delimited accession
 args = [None] * 7
 args[0] = "\"{}\"".format(_bclPath)
 args[1] = "\"{}\"".format(_manifest.lane)
 args[2] = "\"{}\"".format(_manifest.pipeline)
 args[3] = "\"{}\"".format(_manifest.batch)
 args[4] = "\"{}\"".format(_manifest.study)
 args[5] = "\"{}\"".format(",".join(sampleNames))
 args[6] = "\"{}\"".format(_savePath)
 
 if not executeStep("1-bcl-to-bam",args,_executionPath):
  return (False,"Step 1-1 Failed")
 
 print("Step1-1 complete 1/4")
 
 args = [None] * 3
 args[0] = "\"{}\"".format(_savePath)
 args[1] = "\"{}\"".format(_bclPath)
 args[2] = "\"{}\"".format(taglist)
 
 if not executeStep("2-bam-index-decoder",args,_executionPath):
  return (False,"Step 1-2 Failed")
  
 print("Step1-2 complete 2/4")

 args = [None] * 1
 args[0] = "\"{}\"".format(_savePath)
 
 if not executeStep("3-bam-adapter-find",args,_executionPath):
  return (False,"Step 1-3 Failed")
  
 print("Step1-3 complete 3/4")
  
 args = [None] * 1
 args[0] = "\"{}\"".format(_savePath)
 
 if not executeStep("4-splitter",args,_executionPath):
  return (False,"Step 1-4 Failed")
  
  
 print("Step1 complete")
  
 # remove previous steps
 print("Removing intermediate files")
 
 os.remove(os.path.join(_savePath,"step1-1.bam"))
 os.remove(os.path.join(_savePath,"step1-2.bam"))
 os.remove(os.path.join(_savePath,"step1-3.bam"))
 
 return (True,"success")
 
#------------------------------------------------- 
 
def main():
 if len(sys.argv) < 2:
  exitWithError("Missing manifest")
  
 if len(sys.argv) < 3:
  exitWithError("Missing Step 1 scripts folder")
  
 if len(sys.argv) < 4:
  exitWithError("Missing BCL path")
  
 if len(sys.argv) < 5:
  exitWithError("Missing save path")
  
 parser = ASPISManifestParser()
 man = parser.parse(sys.argv[1])
 
 if man == None:
  exitWithError(parser.lastError)
 
 if man[0].batch == None or man[0].batch == "":
   man[0].batch = os.path.basename(sys.argv[1])[0:-4]
 
 result =runStep1(man[0],sys.argv[2],sys.argv[3],sys.argv[4])
 if not result[0]:
  print("Error: {}".format(result[1]))

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Runs step 1 of the Stage 1 Pipeline (BCL-2-VCF)")
 print("usage: python {} [manifest file] [Step 1 scripts folder] [BCL path] [save path]\n".format(sys.argv[0]))
 main();
 print("Bye!")
