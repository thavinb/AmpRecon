import os
import sys
import json
import shutil
import subprocess
import concurrent.futures
import time

#------------------------------------------------- 
def executeStep(_args):
 args = _args
 result = subprocess.run(" ".join(args),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
 if result.returncode != 0:
  return (False,result.stderr.decode("UTF-8"),result.returncode)
 else:
  return (True, result.stdout.decode("UTF-8").strip())

#------------------------------------------------- 

def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()

#------------------------------------------------- 

def qc(_folder):
 
 if not os.path.isdir(_folder):
  _folder = os.path.join(os.getcwd(),_folder)
  
 if not os.path.isdir(_folder):
  return (False, "Could not find cram folder")
  
 crams = sorted([p for p in os.listdir(_folder) if p.endswith(".cram")])
 print("Sample\tline count")
 
 for c in crams:
  lineCount = executeStep(["samtools","view",os.path.join(_folder,c),"|","wc","-l"])
  r = 0
  if lineCount[0]:
   r = lineCount[1]
  else:
   r = "ERROR:{}".format(lineCount[1])
  
  print("{}\t{}".format(c,r))
 
 
 return (True, "done")
   
#------------------------------------------------- 
 
def main():
 if len(sys.argv) < 2:
  exitWithError("Missing cram folder")
  
 result = qc(sys.argv[1])
 if not result[0]:
  print("Error: {}".format(result[1]))

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Performs a line count on all crams found in the selected folder")
 print("usage: python {} [cram folder]\n".format(sys.argv[0]))
 main();
 print("Bye!")
