import re
import os
import sys
import csv
import json
import datetime
import configparser
from enum import Enum


from aspismanifest import ASPISManifest, ASPISManifestParser
from aspisroma import RomaSession

"""
Converts a samplesheet into a ASPIS manifest

[Data]
Sample_ID, Sample_Name, Description, index, index2, Sample_Project

OR

[Data]
sample_ID, sample_name, index, index2, well

REF will be in either description or after _ in Sample_ID

Other columns might be present but will be ignored such as
I7_Index_ID
I5_Index_ID

ASPIS Manifest
---------------
lims id -- Sample_Name
sims id -- Roma ID
index -- Sample count
ref	-- Mapped either from `Description` or from `Sample_ID` `_` sufix
barcode sequence -- index-index2
well - Either `well` or increment
plate -- either generate or take in as an argument

"""

gManifestVersion = "0.1.3"

REQUIRED_FIELD_LIST = ["SAMPLE_ID", "SAMPLE_NAME", "INDEX", "INDEX2"]
REQUIRED_FIELDS = set(REQUIRED_FIELD_LIST)

OPTIONAL_FIELDS = [ "DESCRIPTION", "WELL", "SAMPLE_PROJECT"]

GENERATED_REF = "REF"
GENERATED_INDEX = "POS"
GENERATED_BARCODE = "BARCODE"
GENERATED_ROMA = "ROMA"
GENERATED_PLATE = "PLATE"

PREFERED_REFERENCE_TARGET = OPTIONAL_FIELDS[0]
ALTERNATIVE_REFERENCE_TARGET = REQUIRED_FIELD_LIST[0]

BARCODE_1 = REQUIRED_FIELD_LIST[2]
BARCODE_2 = REQUIRED_FIELD_LIST[3]
BARCODE_FORMAT = "{}-{}"

OPTIONAL_WELL = OPTIONAL_FIELDS[1]
LETTERS = ["A","B","C","D","E","F","G","H"]
WELL_INDICES = [letter+ str(number) for letter in LETTERS for number in range(1,13,1)]
print(len(WELL_INDICES))

SAMPLE_NAME = REQUIRED_FIELD_LIST[1]

REF_LOOKUP = {
    "PFA_GRC1_v1.0":["GRC1"],
    "PFA_Spec":["SP"],
    "PFA_GRC2_v1.0":["GRC2"]
}

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

    args = {"draw":1,"columns":[{"data":"sample_name","name":"","searchable":True,"orderable":True,"search":{"value":"","regex":False}},{"data":"external_id","name":"","searchable":True,"orderable":True,"search":{"value":"","regex":False}}],"order":[{"column":0,"dir":"asc"}],"start":0,"length":25,"search":{"value":_sims,"regex":False}}
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
def samplesheet2manifest(samplesheet_path, save_path):
    print("Checking samplesheet")

    #Load sample sheet
    samplesheet = []
    headers = []
    headers_upper_set = set()


    with open(samplesheet_path) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')

        data_found = False
        for row in reader:
            if data_found:
                if len(headers) == 0:
                    headers = [r.upper() for r in row if len(r.strip()) > 0]
                    headers_upper_set = {r.upper() for r in row if len(r.strip()) > 0}
                else:
                    samplesheet.append({headers[x]: row[x] for x in range(0,len(row)) if len(headers) > x})
                    #print(f"Added {samplesheet[-1]}")
            else:
                if row[0] == "[Data]":
                    data_found = True

    # Validate headers, set mode
    if len(headers_upper_set & REQUIRED_FIELDS) < len(REQUIRED_FIELDS):
        print("Missing required fields: ")
        print(REQUIRED_FIELDS - headers_upper_set)
        return (False, "Missing fields",REQUIRED_FIELDS - headers_upper_set)

    if OPTIONAL_WELL not in headers_upper_set:
        print("Missing Well field, guessing wells")

    #Generate columns
    count = 0
    well_count = 0
    failed_rows = 0
    unique_sample_names = {}
    for row in samplesheet:
        if PREFERED_REFERENCE_TARGET in row:
            row[GENERATED_REF] = row[PREFERED_REFERENCE_TARGET]
        elif '_' in row[ALTERNATIVE_REFERENCE_TARGET]:
            row[GENERATED_REF] = row[ALTERNATIVE_REFERENCE_TARGET].split('_')[1]
        else:
            print(f"Missing reference details for {count + 1}:{row[REQUIRED_FIELD_LIST[0]]}")
            failed_rows += 1

        ref_mapped = False
        if GENERATED_REF in row:
            lookup = row[GENERATED_REF]
            for ref,lookups in REF_LOOKUP.items():
                if lookup in lookups:
                    row[GENERATED_REF] = ref
                    ref_mapped = True
                    break
        if not ref_mapped:
            failed_rows += 1
            print(f"Failed to map reference for {count + 1}:{row[REQUIRED_FIELD_LIST[0]]} lookup:{lookup}")

        row[GENERATED_INDEX] = count + 1
        row[GENERATED_PLATE] = "Sasha-test"
        row[GENERATED_BARCODE] = BARCODE_FORMAT.format(
            row[BARCODE_1],
            row[BARCODE_2]
        )

        if OPTIONAL_WELL not in row:
            row[OPTIONAL_WELL] = WELL_INDICES[well_count]

        if row[SAMPLE_NAME] not in unique_sample_names.keys():
            unique_sample_names[row[SAMPLE_NAME]] = [count]
        else:
            unique_sample_names[row[SAMPLE_NAME]].append(count)

        count += 1
        well_count += 1
        if well_count >= len(WELL_INDICES):
            well_count = 0

    if failed_rows > 0:
        return(False, f"Missing reference details for {failed_rows} samples")

    print("Connecting to ROMA")

    #Create ROMA sessions
    usr = "aspis.pipeline@gmail.com"
    pss = "?9'qn*!K7bZ^:^S"

    romaSessions = {}

    romaUrls = ['https://genre.malariagen.net/'] #'https://spotmalaria.malariagen.net/']
    romaEnums = [RomaInstance.GENRE]
              # RomaInstance.SPOT]

    for j in range(0,len(romaUrls)):
        try:
            romaSessions[romaEnums[j] ] = RomaSession(romaUrls[j], usr,pss)
        except Exception:
            print("Warning unable to contact {} ROMA manifest may contain missing values!".format(romaUrls[j]))
            romaSessions[romaEnums[j] ] = None

    #get unique list of sample names, build roma ids
    #populate samplesheet

    print("Collecting ROMA Ids")
    missing_roma_ids = []
    count = 1
    for sample_name, rows in unique_sample_names.items():
        print(f"Searching {count}/{len(unique_sample_names)}")

        found = False
        romaid = ""
        if sample_name.startswith("Negative_") or sample_name.startswith("3D7_P3"):
            found = True
            romaid = "control"
            print(f"Switching {sample_name} to {romaid}")
        else:
            for roma in romaEnums:
                try:
                    results = romaSessions[roma].query(createQueryObject(sample_name))
                except Exception:
                    results = {}

                if "data" in results and len(results["data"]) > 0:
                    romaid = results["data"][0]["sample_name"]
                    found = True
                    break
        if found:
            for row in rows:
                samplesheet[row][GENERATED_ROMA] = romaid
        count += 1

        if not found:
            print(f"Unable to find Roma Id for {sample_name}")
            missing_roma_ids.append(sample_name)

    if len(missing_roma_ids) > 0:
        print("Failed to find Roma Ids for :-")
        missing_roma_output = "\n".join(missing_roma_ids)
        print(missing_roma_output)
        return (False, f"Missing Roma Ids for {len(missing_roma_ids)} samples", missing_roma_ids)


    print("Buildng Manifest")

    man = ASPISManifest()
    man.version = gManifestVersion
    man.workflow = "LibraryQC"
    man.pipeline = "Sasha - ASPIS"
    man.pipelineVersion = "Alpha"
    man.date = datetime.date.today().isoformat()
    man.lane = 1

    for row in samplesheet:

        if not man.addSample(
            row[SAMPLE_NAME],
            row[GENERATED_ROMA],
            row[GENERATED_INDEX],
            row[GENERATED_REF],
            row[GENERATED_BARCODE],
            row[OPTIONAL_WELL],
            row[GENERATED_PLATE]
            ):
            print("Failed to add {}".format(lims))

    man.save(save_path)

    return(True, "Manifest created")



#-------------------------------------------------
def samplesheet2manifestCmd(args):
    print("\n")
    print(args[0])
    print(len(args[0]) * "=")
    print("Creates a ASPIS manifest from a BCL samplesheet")
    print("Interacts with ROMA to access the missing values")
    print("usage: {} [samplesheet] [save path - optional] \n".format(args[0]))

    if len(args) < 2:
        exitWithError("Missing manifest")

    savePath = os.getcwd()
    if len(args) >= 3:
        savePath = args[2]
    else:
        savePath = os.path.join(
            os.getcwd(),
            f"Manifest_{datetime.date.today().isoformat().replace('/','_')}.csv"
            )

    results = samplesheet2manifest(args[1],savePath)
    if not results[0]:
        exitWithError(results[1])
    print("Bye!")

#-------------------------------------------------
if __name__ == "__main__":
    samplesheet2manifestCmd(sys.argv)
