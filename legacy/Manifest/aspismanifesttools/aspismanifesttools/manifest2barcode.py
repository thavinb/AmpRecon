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

#------------------------------------------------- 
def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()

#------------------------------------------------- 
def createTaglist(_savepath,_manifest):
 print("Creating taglist")
 
 savepath = _savepath
 with open(savepath,"w") as taglist:
  taglist.write("barcode_sequence_1\tbarcode_sequence_2\tbarcode_name\tlibrary_name\n")
  
  if len(sys.argv) >= 4:
   studyId = sys.argv[3]
  else:
   studyId = uuid.uuid4()
  
  controlIndex = 1
  for d in _manifest.data:
   name = d["sims"]
   if name.lower() == "control":
    name = "Control{}".format(controlIndex)
    controlIndex += 1
   bcodes = d["barcode"].split("-")
   taglist.write("{}\t{}\t{}\t{}\n".format(bcodes[0],bcodes[1],d["index"],d["lims"]))
 return savepath
 
#------------------------------------------------- 
def manifest2barcode(manifest,savepath):

 parser = ASPISManifestParser()
 man = parser.parse(manifest)
 
 if man[0] == None:
  exitWithError(man[1])
 else:
  print(createTaglist(savepath,man[0]))

#-------------------------------------------------   
def manifest2barcodeCmd(args):
 print("\n")
 print(args[0])
 print(len(args[0]) * "=")
 print("Creates picard barcode list from a manifest file")
 print("usage: {} [manifest file] [save path]\n".format(args[0]))
 
 if len(args) < 2:
  exitWithError("Missing manifest")
  
 if len(args) < 3:
  exitWithError("Missing save path")
 
 manifest2barcode(args[1],args[2]) 
 
 print("Bye!")
 
 
#------------------------------------------------- 
if __name__ == "__main__":
 manifest2barcodeCmd(sys.argv)
