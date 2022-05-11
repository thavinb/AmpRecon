import json
import sys
import requests
import configparser
import os

try:
 from runStep1 import runStep1
 from runStep2 import runStep2
 from runStep3 import runStep3
 import zipfile
 from zipfile import BadZipFile
 import datetime
 import shutil
 import pysftp
 from google.cloud import storage
 from aspismanifest import ASPISManifest, ASPISManifestParser
 from aspisbackend import ASPISBackendState, ASPISBackendAuthenticator, ASPISBackend
 from aspismanifesttools import manifest2manifest
except Exception as ex:
 print("Import error: {}".format(ex))
 sys.exit(1)
 
#------------------------------------------------- 
URL     = "https://malariagen-aspis.appspot.com/"
WORKER  = "{}-{}".format(sys.argv[0],"1")
VERSION = "0.1.0"

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
 
 if os.path.isfile(os.path.join(os.getcwd(),path)):
  c = configparser.ConfigParser()
  try:
   c.read(path)
  except Exception as ex:
   return (False, ex)
  return (True,c)
 else:
  return (False, "Unable to load / parse config")
  
#------------------------------------------------- 
def exitWithError(_log,error):
 log(_log,"Error: {}".format(error))
 log(_log,"Removing lockfile")
 removeLockFile()
 sys.exit(1)

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
  return (False,"Lock found")
  
 config = loadConfig("launch_stage1.config")
 if not config[0]:
  exitWithError(logFile,config[1])
 config = config[1]
 
 token_auth = ASPISBackendAuthenticator()
 authResults = token_auth.populateFromDict(config["authentication"])
 if not authResults[0]:
  exitWithError(logFile,authResults[1])
 backend = ASPISBackend(URL,WORKER,token_auth)
 
 log(logFile,"Checking for jobs")
 
 #get list of jobs from directory
 jobsPath = os.path.join(os.getcwd(),"jobs")
 jobs = [f for f in os.listdir(jobsPath) if f.endswith(".job")]
 
 if len(jobs) == 0:
  exitWithError(logFile,"No jobs found!")
  
 #flag job as active
 batchId = jobs[0][:jobs[0].find(".")]
 os.rename(os.path.join(jobsPath,jobs[0]),os.path.join(jobsPath,"{}.active".format(batchId)))
 
 #attempt to process as normal
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
  backend.failJob(batchId,e)
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
 
 if manifestObj[0] is None:
  backend.failJob(batchId,manifestObj[1])
  os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.failed".format(batchId)))
  exitWithError(logFile,manifestObj[1])
 manifestObj = manifestObj[0]
 
 backend.reportJobStatus(batchId,"Downloading")
 
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
  backend.failJob(batchId,"Failed to extract BCL")
  os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.failed".format(batchId)))
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
 backend.reportJobStatus(batchId,"Starting Stage 1")
 log(logFile,"Step 1")
 try:
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
     backend.failJob(batchId,result[1])
     os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.failed".format(batchId)))
     exitWithError(logFile,result[1])
   else:
    backend.failJob(batchId,result[1])
    os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.failed".format(batchId)))
    exitWithError(logFile,result[1])
  else:
   backend.failJob(batchId,result[1])
   os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.failed".format(batchId)))
   exitWithError(logFile,result[1])
 except Exception as ex:
  #backend.failJob(batchId,result[1])
  os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.failed".format(batchId)))
  exitWithError(logFile,ex)
  
 log(logFile,"Job has successfully been processed") 
 os.rename(os.path.join(jobsPath,"{}.active".format(batchId)),os.path.join(jobsPath,"{}.done".format(batchId)))
 backend.reportJobStatus(batchId,"Stage-1 complete")
  
 #Clean up step1-4, and step2-12
 log(logFile,"Cleaning up")
 removeFolders = ["step1-4","step2-12"]
 for rFolder in removeFolders:
  rfPath = os.path.join(saveFolder,rFolder)
  if os.path.isdir(rfPath):
   shutil.rmtree(rfPath, ignore_errors = True)

 #Create stage 2 manifest
 log(logFile,"Creating stage 2 manifest")
 stage2ManifestPath = manifest2manifest(manifestObj,os.path.join(saveFolder,"{}_SampleKey.txt".format(batchId)))

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
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Job folder for new jobs and runs stage 1 pipeline if found")
 print("usage: python {}\n".format(sys.argv[0]))
 pollJobs();
 print("\n\nBye!")
