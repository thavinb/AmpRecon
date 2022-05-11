import json
import sys
import requests
import configparser
import os
import datetime
import shutil
import subprocess

from runStep1 import runStep1
from runStep2 import runStep2
from runStep3 import runStep3

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"Manifest","Tools"))
from manifestParser import ASPISManifestParser
from manifestParser import ASPISManifest
from manifest2manifest import createManifest

#-------------------------------------------------
def loadConfig(path = "con.fig"):
 c = configparser.ConfigParser()
 c.read(path)
 return c
 
#------------------------------------------------- 
def activateCondaEnviroment(_name):
 args = ["conda","activate",_name]
 print("Starting conda enviroment {}".format(_name))
 result = subprocess.run(" ".join(args),stderr=subprocess.PIPE,shell=True)
 if result.returncode != 0:
  return (False,"Process failed with error {}\n{}".format(result.returncode,result.stderr.decode("UTF-8")))
 return (True,"Success")
#------------------------------------------------- 
def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()

#------------------------------------------------- 
def runStage1(_batchId,_config = "stage1Poller.config"):
 
 results = activateCondaEnviroment("stage1-pipeline-tools")
 
 # load path from config
 config = loadConfig(_config)
 
 storageConfig = config["storage"]
 savePath = storageConfig["savePath"]
 
 #Check file locations
 bclFolder = os.path.join(config["storage"]["bclpath"],_batchId)
 if not os.path.isdir(bclFolder):
  return (False, "Could not find BCLS for batch {}. Looking for {}".format(_batchId,bclFolder))
  
 manifestPath = os.path.join(config["storage"]["manifestpath"],"{}_manifest.csv".format(_batchId))
 if not os.path.isfile(manifestPath):
  return (False, "Could not find manifest for batch {}. Looking for {}".format(_batchId,manifestPath))
  
 #parse manifest
 print("Loading Manifest")
 parser = ASPISManifestParser()
 manifestObj = parser.parse(manifestPath)
 
 #work around for batch being part of the manifest and study being removed
 manifestObj.batch = _batchId
 manifestObj.study = "Pipeline: {}-{}; Batch: {};".format(manifestObj.pipeline,manifestObj.pipelineVersion,_batchId)
 
 buildSteps = True
 
 #make vcf save path
 print("Creating savepath")
 saveFolder = os.path.join(savePath,_batchId)
 
 if buildSteps:
  if os.path.isdir(saveFolder):
   shutil.rmtree(saveFolder, ignore_errors = True)
  os.mkdir(saveFolder)
 
  #Run tools
  print("Running Stage 1 pipeline")
  print("Step 1")
  result = runStep1(manifestObj,config["scripts"]["step1"],bclFolder,saveFolder)
  if result[0]:
   print("Step 2")
   result = runStep2(manifestObj,config["scripts"]["step2"],config["scripts"]["references"],saveFolder)
   if result[0]:
    print("Step 3")
    result = runStep3(manifestObj,config["scripts"]["step3"],config["scripts"]["references"],saveFolder)
    if result[0]:
     print("Stage 1 pipeline - success!")
    else:
     return (False,result[1])
   else:
    return (False,result[1])
  else:
   return (False,result[1])
  
 stage2Root = config["transfer"]["root"]
 stage2ManPath = os.path.join(stage2Root,"Manifests")
 stage2VCFPath = os.path.join(stage2Root,"Data",_batchId)
  
 #Creating stage 2 manifest
 print("Creating stage 2 manifest")
 stage2ManifestPath = createManifest(stage2ManPath,manifestObj)
 
 #Copy data to stage 2
 print("Copying VCFs to stage 2")
 shutil.copytree(os.path.join(saveFolder,"step3-4"),stage2VCFPath)
 
 return (True,"success")
  
#------------------------------------------------- 
URL     = "https://malariagen-aspis.appspot.com/"
WORKER  = "{}-{}".format(sys.argv[0],"1")
VERSION = "0.1.0"

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Runs ASPIS Stage 1 BCL-TO-VCF pipeline")
 print("usage: python {} [batch id]\n".format(sys.argv[0]))
 
 if len(sys.argv) < 2:
  exitWithError("Missing Batch/Run Id")

 result = runStage1(sys.argv[1])
 
 if not result[0]:
  exitWithError(result[1])

 print("\n\nBye!")
