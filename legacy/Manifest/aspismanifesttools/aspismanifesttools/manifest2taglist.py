import re
import os
import sys
import json
import configparser
import uuid

from aspismanifest import ASPISManifest, ASPISManifestParser

"""
barcode_sequence	barcode_name	library_name	sample_name	description

create a tab delimited file containing the above values

barcode_sequence	- tag1-tag2 (Barcode Sequence)
barcode_name	- Index
library_name	- A number 22906154 !!! - maps to legacy library id. Bambi (decode) needs this but why?
               This is a sanger generated value, does not map to the outside world
               We can make this up <- maybe add SIMS_ID instead?
               doesn't need to be a number
sample_name	- ERS2488229 Accession number <- do we need an accession number?
description - Study Name
"""


def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()
 
#------------------------------------------------- 
def manifest2taglist(manifest, savePath, studyId):
 
 
 if not isinstance(manifest,ASPISManifest):
  parser = ASPISManifestParser()
  man = parser.parse(manifest)
 
  if man[0] == None:
   return(False, man[1])
  man = man[0]
 else:
  man = manifest
  
 print("Creating taglist")
 savepath = os.path.join(savePath,"{}.taglist".format(man.batch))
 with open(savepath,"w") as taglist:
  taglist.write("barcode_sequence\tbarcode_name\tlibrary_name\tsample_name\tdescription\n")
  
  controlIndex = 1
  for d in man.data:
   name = d["sims"]
   if name.lower() == "control":
    name = "Control-{}".format(controlIndex)
    controlIndex += 1
   taglist.write("{}\t{}\t{}\t{}\t{}\n".format(d["barcode"],d["index"],d["lims"],name,studyId))
 return savepath

#------------------------------------------------- 
def manifest2taglistCmd(args):
 
 if len(args) < 2:
  exitWithError("Missing manifest")
  
 if len(args) >= 3:
  savePath = args[2]
 else:
  savePath = os.getcwd()
  
 if len(args) >= 4:
  studyId = args[3]
 else:
  studyId = uuid.uuid4()
  
 print("\n")
 print(args[0])
 print(len(args[0]) * "=")
 print("Creates a taglist from a manifest file")
 print("usage: {} [manifest file] [save path] [study id - optional]\n".format(args[0]))
 print(manifest2taglist(args[1],savePath,studyId))
 
 print("Bye!")
 
#------------------------------------------------- 

if __name__ == "__main__":
 manifest2taglistCmd(sys.argv)
