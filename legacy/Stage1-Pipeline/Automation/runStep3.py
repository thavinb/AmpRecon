import os
import sys
import json
import shutil
import subprocess
import concurrent.futures
import time

from aspismanifest import ASPISManifest, ASPISManifestParser

#------------------------------------------------- 
def executeStep(_target,_args,_execPath):
 args = ["bash",_target +".sh"]
 args.extend(_args)
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

def runStep3ForASample(_executionPath,_refPath, _savePath, _index,_name,_batch):
  
 calls = [ "1-mpilrup",
           "2-call",
           "3-filter1",
           "4-filter2"]
           
 print("Started sample {}-{}".format(_name,_index))
 print("Reference is {}".format(_refPath))

 count = 1
 for c in calls:
  print("Running {} {}/{}".format(c,count,len(calls)))
  
  args = [None] * 3
  args[0] = "\"{}\"".format(_savePath)
  args[1] = "\"{}\"".format(_batch)
  args[2] = "\"{}\"".format(_index)
  
  if c == "1-mpilrup" or "2-call":
   args.append("\"{}\"".format(_refPath))

  if not executeStep(c,args,_executionPath):
   print("Step 3-{} Failed".format(count))
   return False
  
  count += 1
  
 print("Finished sample {}-{}".format(_name,_index))
 
 return True

#------------------------------------------------- 

def runStep3(_manifest,_executionPath,_refPath, _savePath):
 
 # create save paths
 if not os.path.isdir(_savePath):
  return (False,"Save path does not exist")
 
 print("Creating save path")
 stepPaths = ["step3-{}".format(l) for l in range(1,5)]
 for step in stepPaths:
  p = os.path.join(_savePath,step)
  if not os.path.exists(p):
   os.makedirs(p)
 
 #Setup thread pool 
 #lims, sims, accession, index, ref, barcode, library
 sampleTasks = []
 failedTasks = 0
 with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
  for sample in _manifest.data:
   
   refPath = ""
   if sample["ref"] == "PFA_GRC1_v1.0":
    refPath = os.path.join(_refPath,"grc1","Pf_GRC1v1.0")
   elif sample["ref"] == "PFA_Spec":
    refPath = os.path.join(_refPath,"spec","Spec_v1.0")
   elif sample["ref"] == "PFA_GRC2_v1.0":
    refPath = os.path.join(_refPath,"grc2","Pf_GRC2v1.0")
   else:
    print("***********\nNo reference mapping found!\n***********")
   
   sampleTasks.append(executor.submit(runStep3ForASample, _executionPath,refPath, _savePath, sample["index"], sample["sims"],_manifest.batch))
  
  for st in concurrent.futures.as_completed(sampleTasks):
   if not st.result():
    failedTasks += 1
    
  #remove mid process steps
 print("Cleaning up")
 for step in stepPaths[:-1]:
  p = os.path.join(_savePath,step)
  if os.path.exists(p):
   shutil.rmtree(p, ignore_errors = True)
   
 if failedTasks > 0:
  return (False,"{}/{} failed to process!".format(failedTasks,len(sampleTasks)))
 else:
  print("All samples processed!")
  return (True,"success")
   
#------------------------------------------------- 
 
def main():
 if len(sys.argv) < 2:
  exitWithError("Missing manifest")
  
 if len(sys.argv) < 3:
  exitWithError("Missing Step 3 scripts folder")
  
 if len(sys.argv) < 4:
  exitWithError("Missing References path")
  
 if len(sys.argv) < 5:
  exitWithError("Missing save path")
  
 parser = ASPISManifestParser()
 man = parser.parse(sys.argv[1])
 
 if man == None:
  exitWithError(parser.lastError)
 
 result = runStep3(man,sys.argv[2],sys.argv[3],sys.argv[4])
 if not result[0]:
  print("Error: {}".format(result[1]))

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Runs step 3 of the Stage 1 Pipeline (Perform science!)")
 print("usage: python {} [manifest file] [Step 2 scripts folder] [references] [save path]\n".format(sys.argv[0]))
 main();
 print("Bye!")
