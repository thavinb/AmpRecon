from dbUtils import dbAccess
import re
import os
import sys
import json
import configparser
import datetime


def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()
 
 #------------------------------------------------- 

def main():
 
 if len(sys.argv) < 2:
  exitWithError("Missing Run Id")
  
 configPath = "con.fig"
  
 db = dbAccess()
 db.loadConfig(configPath)
 
 sql = """
 SELECT 
 sample.name,sample.supplier_name, qc_result.value,qc_result.assay, 
 qc_result.id_pool_lims, qc_result.units, qc_result.labware_purpose, 
 qc_result.qc_type,qc_result.last_updated
 FROM sample left join qc_result on qc_result.id_sample_tmp = sample.id_sample_tmp
 WHERE
 sample.name in (select 
				sample.name from sample, iseq_flowcell , iseq_product_metrics
				where iseq_flowcell.id_sample_tmp = sample.id_sample_tmp 
				AND iseq_product_metrics.id_iseq_flowcell_tmp = iseq_flowcell.id_iseq_flowcell_tmp 
				AND iseq_product_metrics.id_run = %s
    AND sample.supplier_name != "control")
 and qc_result.labware_purpose = 'SWGA Stock' and (qc_result.assay = 'QuantEssential v0.1' or qc_result.assay = 'Volume Check 1.0')
 order by qc_result.last_updated DESC
 """
 
 results = db.query(sql,(sys.argv[1],))
 
 if results and len(results) > 0:
  
  sample = {}
  
  for r in results:
   if r[0] not in sample:
    if r[3] == "QuantEssential v0.1":
     sample[r[0]] = [r[1],r[4],str(r[8]),"{} {}".format(r[2],r[5]),"-"]
    elif r[3] == "Volume Check 1.0":
     sample[r[0]] = [r[1],r[4],str(r[8]),"-","{} {}".format(r[2],r[5])]
    else:
     print("Unknown assay type found {}!".format(r[3]))
   else:
    if r[3] == "QuantEssential v0.1" and sample[r[0]][3] == "-":
     sample[r[0]][3] = "{} {}".format(r[2],r[5])
    elif r[3] == "Volume Check 1.0" and sample[r[0]][4] == "-":
     sample[r[0]][4] = "{} {}".format(r[2],r[5])
    else:
     print("Failed to set extra assay type! {}: [3]={} [4]={}".format(r[0],sample[r[0]][3],sample[r[0]][3]))
  
  with open("{}_sampleVolumes.csv".format(sys.argv[1]),"w") as svfile:
   svfile.write(", ".join(["Name","Supplier","Pool","Last Update","Concentration","Volume"]) + "\n")
   for k,v in sample.items():
    svfile.write("{}, {}\n".format(k,", ".join(v)))
  
  
 else:
  exitWithError("Could not find run!")

 

#------------------------------------------------- 
if __name__ == "__main__":

 print("\n")
 print(sys.argv[0])
 print(len(sys.argv[0]) * "=")
 print("Creates a sample volume csv for a given run")
 print("usage: python {} [run id]\n".format(sys.argv[0]))
 main();
 print("Bye!")
