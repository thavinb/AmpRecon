import re
import os
import sys
import json
import configparser
import datetime

from aspisdb.utils.v0_1_0 import dbAccess
from aspismanifest import ASPISManifest


gManifestVersion = "0.1.3"

def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()


#-------------------------------------------------
def getMLWHSql():
 return """SELECT sample.name,sample.accession_number,iseq_flowcell.primer_panel,
 iseq_flowcell.tag_index,sample.supplier_name,iseq_product_metrics.id_run,
 iseq_flowcell.tag_sequence,iseq_flowcell.tag2_sequence,
 iseq_flowcell.legacy_library_id,
 stock_resource.labware_coordinate, stock_resource.labware_human_barcode
 FROM iseq_product_metrics, iseq_flowcell, sample, stock_resource
 WHERE iseq_product_metrics.id_run = %s AND
 iseq_product_metrics.id_iseq_flowcell_tmp = iseq_flowcell.id_iseq_flowcell_tmp
 and iseq_flowcell.id_sample_tmp = sample.id_sample_tmp
 and stock_resource.id_sample_tmp = sample.id_sample_tmp
 order by iseq_flowcell.tag_index"""


#-------------------------------------------------
def dumpSanger(run,savePath,configpath):
 db = dbAccess()
 db.loadConfig(configpath)

 sql = getMLWHSql()

 results = db.query(sql,(run,))

 if results and len(results) > 0:
  with open(savePath,"w") as dump:
   dump.write("name,supplier_name,tag_index,primer_panel,tag_sequence,tag2_sequence,labware_coordinate, labware_human_barcode\n")
   for r in results:
    lims, accession,ref, index, sims, runid, tag1, tag2,library, well,plate = r
    dump.write("{},{},{},{},{},{},{},{}\n".format(lims,sims,index, ref,tag1, tag2,well,plate))
  return (True,"dumped")
 else:
  return (False,"Failed to access results")

#-------------------------------------------------


def sanger2Manifest(run,pipe,semver,configpath,savePath = ""):

 db = dbAccess()
 db.loadConfig(configpath)

 sql = getMLWHSql()

 results = db.query(sql,(run,))

 if results and len(results) > 0:

  #Validate
  invalid = []
  romaRegEx    = re.compile("(SPT|RCN|VVX|VBS)\d+")
  controlRegEx = re.compile("([cC]ontrol|(CAS|CHP|CNG|CLP|ING|ILP)\d+)")

  for r in results:
   romaId = r[4]
   if romaId is None or romaRegEx.fullmatch(romaId) is None:
    if romaId is None or controlRegEx.fullmatch(romaId) is None:
     invalid.append(romaId)


  if len(invalid) > 0:
   return (False,"{} invalid supplier Ids found in MLWH".format(len(invalid)),"bad-roma")

  #create file
  manifestName = "{}_manifest.csv".format(run)

  if len(savePath) > 0:
   manifestName = savePath

  man = ASPISManifest()
  man.version = gManifestVersion
  man.workflow = "LibraryQC"
  man.pipeline = pipe
  man.pipelineVersion = semver
  man.date = datetime.date.today().isoformat()
  man.lane = 1

  for r in results:
   lims, accession,ref, index, sims, runid, tag1, tag2,library, well,plate = r
   if sims == "Control":
    sims = sims.lower()

   if not man.addSample(lims,sims,index, ref,"{}-{}".format(tag1, tag2),well,plate):
    print("Failed to add {}".format(lims))

  results = man.save(manifestName)

  if results[0]:
   return (True, "Manifest created: {}".format(manifestName))
  else:
   return (False, "Failed to create manifest - {}".format(results[1]))
 else:
  return (False,"Could not find run!")


#-------------------------------------------------
def sanger2ManifestCmd(args):
 print("\n")
 print(args[0])
 print(len(args[0]) * "=")
 print("Create a Stage-1 Pipeline manifest file")
 print("usage: {} [run id] [pipeline] [pipeline version] [config path]\n".format(args[0]))

 if len(args) < 2:
  exitWithError("Missing Run Id")

 if len(args) < 3:
  exitWithError("Missing Pipeline")

 if len(args) < 4:
  exitWithError("Missing Pipeline version")

 config = "~/mlwh-config"
 if len(args) >= 5:
  config = args[4]

 results = sanger2Manifest(args[1],args[2],args[3],config);

 if results[0]:
  print("Success: {}".format(results[1]))
 else:
  exitWithError(results[1])
 print("Bye!")

#-------------------------------------------------
if __name__ == "__main__":
 sanger2ManifestCmd(sys.argv)
