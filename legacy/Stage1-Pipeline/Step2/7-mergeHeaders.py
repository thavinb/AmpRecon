import sys
import time
import os
import shutil
import subprocess
import re
from signal import signal, SIGPIPE, SIG_DFL

"""
where SN field in SQ header record matches one in the given dict file,
replace that SQ record with that in the dict file,
but propagate any 'AH' fields provided by the aligner"
"""


def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit(2)
 
#------------------------------------------------- 
 
def querySamtools(target):

 result = subprocess.run(["samtools","view","-H",target],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 if result.returncode == 0:
  return result.stdout.decode("UTF-8")
 else:
  print("Error: {}".format(result.stderr.decode("UTF-8")))
  return None
  
#------------------------------------------------- 

def getHeaderValue(line,target):
 if target in line:
  start = line.index(target)
  end   = line.index("\t",start)
  return line[start:end]
 return ""
#------------------------------------------------- 

def main():
 if len(sys.argv) < 2:
  exitWithError("Index id not found")
 
 if len(sys.argv) < 3:
  exitWithError("Missing dict")
  
 if len(sys.argv) < 4:
  exitWithError("Missing save path")
  
 if len(sys.argv) < 5:
  exitWithError("Missing batch number")
  
 signal(SIGPIPE, SIG_DFL) 
 regex = re.compile("@SQ.*\tSN:([^\t]+)")

 old = querySamtools("{}/step2-3/{}_1#{}.cram".format(sys.argv[3],sys.argv[4],sys.argv[1]))
 new = querySamtools("{}/step2-6/{}_1#{}.bam".format(sys.argv[3],sys.argv[4],sys.argv[1]))
 
 lastPGID = ""
 filteredOldArray = []
 for l in old.split("\n"):
  
  filteredOldArray.append(l)
  if l.startswith("@PG") and "ID:" in l:
   lastPGID = getHeaderValue(l,"ID:").replace("ID","\tPP")

 old = "\n".join(filteredOldArray) 
 
 lookup = {}
 with open("{}.dict".format(sys.argv[2]),"r") as fd:
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
 

#------------------------------------------------- 
if __name__ == "__main__":
 main();
else:
 print(__name__)
