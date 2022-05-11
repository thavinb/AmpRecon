import re
import os
import sys
import json
import configparser
import uuid

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

def formatBarcode(_bar):
 if len(_bar) == 16:
  return "{}-{}".format(_bar[0:8],_bar[8:])
 else:
  return _bar


#------------------------------------------------- 

def scanBarcodes(_path):
 if not os.path.isdir(_path):
  exitWithError("{} is not a directory".format(_path))
  
 results = [p for p in os.listdir(_path) if p.endswith(".txt")]
 
 samples = {}
 samplePool = []
 
 for res in results:
  with open(os.path.join(_path,res),"r") as rin:
   for l in rin:
    columns = l.split('\t')
    if len(columns) != 5:
     print("Skipping {}, unexpected format".format(res))
     break
    columns[2] = columns[2].upper()
    if columns[1].upper() == "Y":
     if columns[2] not in samples:
      samples[columns[2]] = 1
      if columns[2] in samplePool:
       samplePool.remove(columns[2])
     else:
      samples[columns[2]] += 1
    else:
     if("." not in columns[2] 
     and len(columns[2]) > 0 
     and columns[2] not in samplePool 
     and columns[2] not in samples):
      samplePool.append(columns[2])
      
 if len(samples) == 0:
  print("No barcodes found!")
 else:
  filteredOut = 0
  filterPoint = 5
  for k,v in samples.items():
   if v > filterPoint:
    print("{} found {} times".format(formatBarcode(k),v))
   else:
    filteredOut += 1
  if filteredOut > 0:
   print("{} barcodes where found less than {} times and have not been printed".format(filteredOut,filterPoint)) 
    
  print("{} barcodes found".format(len(samples)))
 
 if len(samplePool) > 0:
  print("\nMissing barcodes")
  for s in samplePool:
   print(formatBarcode(s))
  print("{} Samples barcodes".format(len(samplePool)))
 
#------------------------------------------------- 

def barcodeReviewCMD():
 
 if len(sys.argv) < 2:
  exitWithError("Missing barcode path")

 scanBarcodes(sys.argv[1])
 
#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Scans the output of Picard-ExtractIlluminaBarcodes and outputs results")
 print("usage: python {} [barcode path]\n".format(sys.argv[0]))
 barcodeReviewCMD();
 print("Bye!")
