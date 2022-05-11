import json
import sys
import requests
import configparser
import os
import zipfile
from zipfile import BadZipFile
import datetime
import shutil

import pysftp

from google.cloud import storage

from runStep1 import runStep1
from runStep2 import runStep2
from runStep3 import runStep3

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"Manifest","Tools"))
from manifestParser import ASPISManifestParser
from manifestParser import ASPISManifest
from manifest2manifest import createManifest

#-------------------------------------------------
def createLockFile():
 """
 Creates a file in the CWD to show that it is currently executing
 """
 lock = open(datetime.datetime.today().strftime("%Y-%m-%d-%H-%M-%S.lock"),"w")

#-------------------------------------------------
def removeLockFile():
 """
 Removes the file created be createLockFile
 """
 lockFiles = [f for f in os.listdir(os.getcwd()) if os.path.isfile(os.path.join(os.getcwd(),f)) and f.endswith(".lock")]
 for l in lockFiles:
  os.remove(l)
  
#-------------------------------------------------
def hasLockFile():
 """
 Checks to see if a lock file exists
 """
 for f in os.listdir(os.getcwd()):
  if f.endswith(".lock"):
   return True
 return False

#-------------------------------------------------
def loadConfig(path = "con.fig"):
 c = configparser.ConfigParser()
 c.read(path)
 return c
 
#------------------------------------------------- 
def exitWithError(_log,error):
 log(_log,"Error: {}".format(error))
 log(_log,"Removing lockfile")
 removeLockFile()
 sys.exit()

#------------------------------------------------- 
def auth_token(_config):
 
 backendConfig = _config["backend"] 
 token_url = backendConfig['token_url']
 
 args = {}
 for key in backendConfig:
  args[key] = backendConfig[key]
 
 token_response = requests.get(token_url, args)
 if token_response.status_code == 401:
     return None
 token = token_response.text.split('=')
 token = token[1].split('&')[0]
 return token
 
#------------------------------------------------- 
def unzipBCL(_bcl,_batch,_savepath):
 with zipfile.ZipFile(_bcl,"r") as job: 
  if job.testzip() == None:
   bclFolder = os.path.join(_savepath,_batch)
   if os.path.isdir(bclFolder):
    shutil.rmtree(bclFolder, ignore_errors = True)
   os.mkdir(bclFolder)
       
   job.extractall(bclFolder)
   return bclFolder
 return None
#------------------------------------------------- 
def log(_log,_msg):
 _msg = "{} - {}".format(datetime.datetime.today().strftime("%H:%M:%S"),_msg)
 _log.write(_msg + "\n")
 print(_msg)

#-------------------------------------------------

def createLogFile():
 if not os.path.isdir("logs"):
  os.mkdir("logs")
 
 logFile = open(os.path.join("logs",datetime.datetime.today().strftime("%Y-%m-%d-%H-%M-%S.log")),"w")
 log(logFile,"Running {} {}\n".format(sys.argv[0],VERSION))
 return logFile

#------------------------------------------------- 
def findBclRoot(_folder):
 files = [f for f in os.listdir(_folder) if os.path.isfile(os.path.join(_folder,f))]
 if "RunInfo.xml" in files:
  return _folder
 else:
  folders = [f for f in os.listdir(_folder) if os.path.isdir(os.path.join(_folder,f))]
  if len(folders) == 0:
   return None
   
  for f in folders:
   result = findBclRoot(os.path.join(_folder,f))
   if result is not None:
    return result
 
#------------------------------------------------- 

def pollJobs():
 
 logFile = createLogFile()
 
 log(logFile,"Checking for lock")
 if hasLockFile():
  log(logFile,"Lock found, exiting")
  return False

 config = loadConfig("stage1Poller.config")
 
 token = auth_token(config)
 
 if token is None:
  exitWithError(logFile,"Failed to authenticate")
 
 log(logFile,"Checking for jobs")
 jsJob = json.loads(requests.post(URL+'work?id=' + WORKER, headers={'Authorization': 'Bearer '+token}).text)
 if "batch_id" not in jsJob:
  exitWithError(logFile,"No jobs found!")

 batchId = jsJob["batch_id"]

 log(logFile,"Accessing Batch {}".format(batchId))

 bclPath      = "{}/BCL".format(batchId)
 manifestPath = "{}/manifest".format(batchId)

 storageConfig = config["storage"]
 gcpBucket = storageConfig["bucket"]
 manifestSavePath = storageConfig["manifestPath"]
 bclSavePath = storageConfig["bclPath"]
 savePath = storageConfig["savePath"]

 sc = storage.Client()
 try:
  bucket = sc.get_bucket(gcpBucket)
 except Exception as e:
  exitWithError(logFile,e)
 
 #create lock file
 log(logFile,"Creating lock file")
 createLockFile()
 
 #Save manifest
 blob = bucket.get_blob(manifestPath)
 manifestSave = os.path.join(manifestSavePath,"{}.csv".format(batchId))
 log(logFile,"Saving manifest to {}".format(manifestSave))
 with open(manifestSave,"wb") as file_obj:
  blob.download_to_file(file_obj)
  
 #parse manifest
 log(logFile,"Loading Manifest")
 parser = ASPISManifestParser()
 manifestObj = parser.parse(manifestSave)
  
 #Save BCL
 blob = bucket.get_blob(bclPath)
 bclSave = os.path.join(bclSavePath,"{}.zip".format(batchId))
 log(logFile,"Saving bcl to {}".format(bclSave))
 with open(bclSave,"wb") as file_obj:
  blob.download_to_file(file_obj)
  
 #Unzip BCL
 log(logFile,"Unpacking BCL")
 #bclFolder = os.path.join(bclSavePath,batchId)
 bclFolder = unzipBCL(bclSave,batchId,bclSavePath)
 if bclFolder is None:
  exitWithError(logFile,"Failed to extract BCL")
  
 #Remove BCL zip
 log(logFile,"Removing BCL zip")
 os.remove(bclSave)
 
 #find root folder
 log(logFile,"Finding BCL root")
 bclFolder = findBclRoot(bclFolder)
 log(logFile, "BCL root is {}".format(bclFolder))
 
 #make vcf save path
 log(logFile,"Creating savepath")
 saveFolder = os.path.join(savePath,batchId)
 if os.path.isdir(saveFolder):
  shutil.rmtree(saveFolder, ignore_errors = True)
 os.mkdir(saveFolder)
 
 #work around for batch being part of the manifest and study being removed
 manifestObj.batch = batchId
 manifestObj.study = "Pipeline: {}-{}; Batch: {};".format(manifestObj.pipeline,manifestObj.pipelineVersion,batchId)

 #Run tools
 #log(logFile,"Skipping Stage 1 pipeline")
 log(logFile,"Running Stage 1 pipeline")
 log(logFile,"Step 1")
 result = runStep1(manifestObj,config["scripts"]["step1"],bclFolder,saveFolder)
 if result[0]:
  log(logFile,"Step 2")
  result = runStep2(manifestObj,config["scripts"]["step2"],config["scripts"]["references"],saveFolder)
  if result[0]:
   log(logFile,"Step 3")
   result = runStep3(manifestObj,config["scripts"]["step3"],config["scripts"]["references"],saveFolder)
   if result[0]:
    log(logFile,"Stage 1 pipeline - success!")
   else:
    exitWithError(logFile,result[1])
  else:
   exitWithError(logFile,result[1])
 else:
  exitWithError(logFile,result[1])
  
 log(logFile,"Job has successfully been processed") 
  
 #Clean up step1-4, and step2-12
 log(logFile,"Cleaning up")
 removeFolders = ["step1-4","step2-12"]
 for rFolder in removeFolders:
  rfPath = os.path.join(saveFolder,rFolder)
  if os.path.isdir(rfPath):
   shutil.rmtree(rfPath, ignore_errors = True)

 #Create stage 2 manifest
 log(logFile,"Creating stage 2 manifest")
 stage2ManifestPath = createManifest(saveFolder,manifestObj)
 
 #Upload to stage 2
 transferConfig = config["transfer"]
 with pysftp.Connection(transferConfig["host"], username=transferConfig["user"], private_key=transferConfig["pam"]) as sftp:

  #Upload VCFs
  log(logFile,"Upload VCFs")
  vcfPath = os.path.join(transferConfig["root"] ,"Data",batchId)
  stepPath = os.path.join(saveFolder,"step3-4")
  
  #Create folder
  if sftp.isdir(vcfPath):
   sftp.rmdir(vcfPath)
   
  sftp.mkdir(vcfPath)
  sftp.put_d(stepPath,vcfPath) 
  
  log(logFile,"Upload stage 2 manifest")
  with sftp.cd(os.path.join(transferConfig["root"] ,'Manifests')):
   sftp.put(stage2ManifestPath) 
 
 #Unlock
 log(logFile,"Removing lockfile")
 removeLockFile()
 
 log(logFile,"Done")
 
  
#------------------------------------------------- 
URL     = "https://malariagen-aspis.appspot.com/"
WORKER  = "{}-{}".format(sys.argv[0],"1")
VERSION = "0.1.0"

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Polls ASPIS backend for new jobs and runs stage 1 pipeline if found")
 print("usage: python {}\n".format(sys.argv[0]))
 pollJobs();
 print("\n\nBye!")
