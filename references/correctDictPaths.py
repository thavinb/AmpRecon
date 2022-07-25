import os
import sys
import json

#------------------------------------------------- 

def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()

#------------------------------------------------- 

def updateDict(_path):
 
 if os.path.isdir(_path):
  path = _path
 else:
  path = os.path.join(os.getcwd(),_path)
 
 files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f)) ]

 fastaPath = ""
 #find fasta
 for f in files:
  if f.lower().endswith(".fasta"):
   fastaPath = os.path.join(path,f)
   break
 
 dictPath = ""
 #find dict
 for f in files:
  if f.lower().endswith(".dict"):
   dictPath = os.path.join(path,f)
   break
   
 if len(fastaPath) == 0:
  print("Failed to find fasta")
  return False
  
 if len(dictPath) == 0:
  print("Failed to find dict")
  return False
 
 #scan dict
 dictLines = []
 with open(dictPath,"r") as dictIn:
  for line in dictIn:
   dictLines.append(line)
 
 target = "UR:file:"
 oldFilePath = ""
 #find current location
 for l in dictLines:
  if target in l:
   columns = l.split('\t')
   for c in columns:
    if target in c:
     oldFilePath = c[len(target):]
     break
   break
 
 if len(oldFilePath) == 0:
  return False
 
 #replace file path
 for i in range(0,len(dictLines)):
  dictLines[i] = dictLines[i].replace(oldFilePath,fastaPath)
 
 #rewrite dicts
 with open(dictPath,"w") as dictOut:
  for line in dictLines:
   dictOut.write(line)
 
 return True
#------------------------------------------------- 
 
def main():
 #list folders
 folders = [f for f in os.listdir() if os.path.isdir(f)]
 #create
 for f in folders:
  print("Processing {}".format(f))
  if updateDict(f):
   print("Directory updated")
  else:
   print("Failed to update")

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Corrects the dictionary fasta file paths to point to the new fasta locations")
 print("Expects to be ran in the root folder containing grc1,grc2, and spec folders")
 print("usage: python {}\n".format(sys.argv[0]))
 main();
 print("Bye!")
