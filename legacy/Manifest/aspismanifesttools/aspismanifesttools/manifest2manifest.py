import re
import os
import sys
import json
import configparser
from enum import Enum

from aspismanifest import ASPISManifest, ASPISManifestParser
from aspisroma import RomaSession

"""
Stage 2 manifest format - tab deliimited
Sanger_ID	Index	Plate_ID	Well	Supplier_ID	External_ID	Date	Location	Country	Study
M         M     M        M    M           R           R    R        R       R

"""
#-------------------------------------------------

class RomaInstance(Enum):
 SPOT  = 1
 GENRE = 2
 VIVAX = 3
 VOBS  = 4

#-------------------------------------------------
def MapSIMSToURL(_sims):
 if _sims.startswith("SPT"):
  return RomaInstance.SPOT
 elif _sims.startswith("RCN"):
  return RomaInstance.GENRE
 elif _sims.startswith("VVX"):
  return RomaInstance.VIVAX
 elif _sims.startswith("VBS"):
   return RomaInstance.VOBS
 else:
  return None

#-------------------------------------------------
def exitWithError(error):
 print("Error: {}".format(error))
 sys.exit()
 
#-------------------------------------------------
def createQueryObject(_sims):
 
 if not _sims.isdigit():
  _sims = _sims[3:]
  
  if _sims.startswith("0"):
   _sims = _sims[1:]
   
 args = {"draw": 1,
          "columns": [
          {"data": "id", "name": "", "searchable": True, "orderable": True,"search": {"value": "", "regex": False}},
          {"data": "collection_date", "name": "", "searchable": True, "orderable": True,"search": {"value": "", "regex": False}}
          ], "order": [{"column": 0, "dir": "asc"}
          ], "start": 0, "length": 1,
          "search": {"value": _sims, "regex": False}}
 return args

#------------------------------------------------- 
def getCollectionValue(_tag,_target):
 if "collection_" not in _target:
  _target = "collection_" + _target
  
 start = _tag.find(_target)
 keyValue = _tag[start:_tag.find(",",start)]
 kvpair = keyValue.split(":")
 
 return kvpair[1].replace("\"","").strip()
 
 
#------------------------------------------------- 
def manifest2manifest(manifest,savepath):
 print("Creating manifest")
 
 if not isinstance(manifest,ASPISManifest):
  parser = ASPISManifestParser()
  man = parser.parse(manifest)
 
  if man[0] == None:
   return(False, man[1])
  manifest = man[0]

 print("Batch set to {}".format(manifest.batch))
 
 sample_columns = ['external_id',  'location_name', 'country_id','study_name']

 if savepath is None:
  savepath = "{}_SampleKey.txt".format(manifest.batch)
  
 with open(savepath,"w") as taglist:
  taglist.write("Sanger_ID\tIndex\tPlate_ID\tWell\tSupplier_ID\tExternal_ID\tDate\tLocation\tCountry\tStudy\n")
  
  #Create ROMA sessions
  usr = "aspis.pipeline@gmail.com"
  pss = "?9'qn*!K7bZ^:^S"
  
  romaSessions = {}
  
  romaUrls = ['https://spotmalaria.malariagen.net/',
              'https://genre.malariagen.net/',
              'https://vivax.malariagen.net/',
              'https://vobs.malariagen.net/']
  romaEnums = [RomaInstance.SPOT,
               RomaInstance.GENRE,
               RomaInstance.VIVAX,
               RomaInstance.VOBS]

  for j in range(0,len(romaUrls)):
   try:
    romaSessions[romaEnums[j] ] = RomaSession(romaUrls[j], usr,pss)
   except Exception:
    print("Warning unable to contact {} ROMA manifest may contain missing values!".format(romaUrls[j]))
    romaSessions[romaEnums[j] ] = None

  vPlate_ID    = None #	- ?
  vWell	       = None #	- ?
  vSupplier_ID	= None # - R - got sample_name
  vExternal_ID = None # - R - got external_id
  vDate	       = None # - R   got tags.collection_[day,month,year] <- optional!
  vLocation    = None # - R - got location.text
  vCountry     = None # - R - got location.country.short_name
  vStudy	      = None # - R - got manifest.study.title
  
  missingValue = "-"
  
  controlIndex = 1
  for d in manifest.data:
   #create url based on sims name
   romaType = MapSIMSToURL(d["sims"])
   
   vPlate_ID    = d["plate"] 
   vWell	       = d["well"]
   vSupplier_ID	= missingValue
   vExternal_ID = missingValue 
   vDate	       = missingValue 
   vLocation    = missingValue 
   vCountry     = missingValue
   vStudy	      = missingValue
   
   if romaType is not None and romaSessions[romaType] is not None:
    
    try:
     results = romaSessions[romaType].query(createQueryObject(d["sims"]))
    except Exception:
     results = {}
    
    if "data" in results and len(results["data"]) > 0:
     r = results["data"][0]
     
     #tags.collection_[day,month,year] <- optional!
     if "collection_date" in r and r["collection_date"] is not None:
      vDate = r["collection_date"]
     elif "tags" in r and r["tags"] is not None:
      if "collection_day" in r["tags"] and \
         "collection_month" in r["tags"] and \
         "collection_year" in r["tags"]:
       
       vDate = "{}/{}/{}".format(getCollectionValue(r["tags"],"year"),
                                 getCollectionValue(r["tags"],"month"),
                                 getCollectionValue(r["tags"],"day"))
 
     if vDate == missingValue:
      print("Unable to get collection date for {}".format(d["sims"]))

     vSupplier_ID	= r["sample_name"]
     vExternal_ID = r["external_id"]
     vLocation    = r["location"]["text"]
     vCountry     = r["location"]["country"]["short_name"]
     vStudy	      = "{}:{}".format(r["manifest"]["study"]["project_code"],r["manifest"]["study"]["title"])
     
    else:
     print("Failed to collect ROMA data for {}".format(d["sims"]))
     print(results)
    
   taglist.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(d["lims"],d["index"],vPlate_ID,vWell,d["sims"],vExternal_ID,vDate,vLocation,vCountry,vStudy))
 return savepath
 
#------------------------------------------------- 
def manifest2manifestCmd(args):
 print("\n")
 print(args[0])
 print(len(args[0]) * "=")
 print("Creates a stage 2 manifest from a stage 1 manifest file")
 print("Interacts with ROMA to access the missing values")
 print("usage: {} [manifest file] [save path - optional] \n".format(args[0]))
 
 if len(args) < 2:
  exitWithError("Missing manifest")
  
 savePath = os.getcwd()
 if len(args) >= 3:
  savePath = args[2]
 
 results = manifest2manifest(args[1],savePath)
 if not results[0]:
  exitWithError(results[1])
 print("Bye!")
 
#------------------------------------------------- 
if __name__ == "__main__":
 manifest2manifestCmd(sys.argv)
