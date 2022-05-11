import re
import os
import sys
import json
import configparser
import datetime

from aspisdb.utils.v0_1_0 import dbAccess

gManifestVersion = "0.1.0"

"""
[Header]
Version,	0.1.0
Workflow,	LibraryQC
Date,	2018-04-23
Platform Unit,	???
Study Name,	???
Lane,	1

[Data]
LIMS ID,	SIMS ID,	Index,	Ref, 	Barcode Sequence

Getting None as a value for Accession number for job 27835
"""

def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()

#------------------------------------------------- 

def injectAccessionIntoGRC(manifest,config):
   
 db = dbAccess()
 db.loadConfig(config)
 
 # load GRC
 if not os.path.isfile(sys.argv[1]):
  exitWithError("Unable to find GRC")
 
 print("Reading GRC") 
 grcLines = []
 with open(sys.argv[1],"r") as grcf:
  for line in grcf:
   grcLines.append(line)
   
 # todo add GRC format validation
   
 # get accessions for run id
 runId = os.path.basename(sys.argv[1])
 runId = runId[0:runId.find("_")]
 
 print("Looking for run {}".format(runId))
 
 sql = "SELECT sample.accession_number,sample.supplier_name "
 sql+= "FROM iseq_product_metrics, sample, iseq_flowcell "
 sql+= "WHERE iseq_product_metrics.id_run = %s AND "
 sql+= "iseq_product_metrics.id_iseq_flowcell_tmp = iseq_flowcell.id_iseq_flowcell_tmp "
 sql+= "and iseq_flowcell.id_sample_tmp = sample.id_sample_tmp "
 sql+= "order by iseq_flowcell.tag_index"
 
 results = db.query(sql,(runId,))
 
 if results is None or len(results) == 0:
  exitWithError("Unable to find run id {}".format(runId))
 
 accessionMap = {}
 
 for s in results:
  if s[0] not in accessionMap:
   accessionMap[s[1]] = s[0]
   
 # Create new GRC
 outPath = sys.argv[1]
 outPath = outPath[:-4] + "-with-accession.txt"
 
 print("Saving to {}".format(outPath))

 with open(outPath,"w") as grcn:
  
  grcn.write("Accession Number\t" + grcLines[0])

  for sample in grcLines[1:]:
   sims = sample[0:sample.find("\t")]
  
   lineFormat = "{}\t{}"
   if sims in accessionMap:
    grcn.write(lineFormat.format(accessionMap[sims],sample))
   else:
    print("ERROR! - Could not find sample {} in MLWH".format(sims))
    grcn.write(lineFormat.format("-",sample))
  
#------------------------------------------------- 
def injectAccessionIntoGRCCMD():
 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Injects Accession nuymbers into a GRC")
 print("usage: python {} [GRC path] [config - optinal]\n".format(sys.argv[0]))
 
 if len(sys.argv) < 2:
  exitWithError("Missing GRC")
  
 configPath = "con.fig"
 if len(sys.argv) >= 3:
  configPath = sys.argv[2]
 
 injectAccessionIntoGRC(sys.argv[1],configPath);
 print("Bye!")

#------------------------------------------------- 
if __name__ == "__main__":
 injectAccessionIntoGRCCMD()

