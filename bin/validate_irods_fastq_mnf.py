#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

from ast import Assert
import pandas as pd
from sys import argv, exit

# load manifest
mnf_flpath = argv[1]
panel_set_flpath = argv[2]
execution_mode = argv[3]
df = pd.read_csv(mnf_flpath, sep='\t')

# Set global requirements depending which execution mode is selected
if execution_mode=='irods':
    REQUIRED_HEADERS = ["sample_id", "primer_panel","irods_path"]
    REQUIRED_FORMAT = ".cram"
    name_path_column = "irods_path"
elif execution_mode=='fastq':
    REQUIRED_HEADERS = ["sample_id", "primer_panel","fastq_path"]
    REQUIRED_FORMAT = tuple(['.fastq','.fq', '.fastq.gz', '.fq.gz'])
    name_path_column = "fastq_path"


pnl_df = pd.read_csv(panel_set_flpath, sep=',')
REQUIRED_PANELS = pnl_df["panel_name"].to_list()
ERRORS_FOUND = 0

# check if required headers are present
print("@ checking column headers...")
columns_lst = df.columns.tolist()
miss_headers = []
for rq_h in REQUIRED_HEADERS:
    found = False
    for col in columns_lst:
        if rq_h == col:
            found=True
            break
    if found == False:
        miss_headers.append(rq_h)
try:
    assert(len(miss_headers) == 0)
    print("    > PASS")

except(AssertionError):
    print("ERROR: The following expected headers were no found: ")
    print(f"      {miss_headers}")
    ERRORS_FOUND +=1

# check if expected patterns for irods_path were provided
# is it only cram files?
print(f"@ checking if only expected format ['{REQUIRED_FORMAT}'] was requested...")
def isItRequiredFormat(row):
    return row.endswith(REQUIRED_FORMAT)

#name_path_column is either irods_path or fastq_path depending on the execution mode
valid_format_idx = df[name_path_column].apply(isItRequiredFormat)
incorrect_format_rows = df.loc[~valid_format_idx]

try:
    assert(len(incorrect_format_rows) == 0)
    print("    > PASS")
except:
    print("ERROR: The manifest contains paths to files with unsupported formats. Only CRAM (execution modes: irods) and FASTQ files (execution mode: fastq) are supported.")
    print(f"{incorrect_format_rows[REQUIRED_HEADERS]}")
    ERRORS_FOUND +=1

# check if WG_LANE_SAMPLE_ID is unique
print("@ check for pipeline internal ids uniqueness...")
def genPipelineInternalId(row, execution_mode):
    if(execution_mode=='irods'):
        simple_name = row[name_path_column].split("/")[-1].split(".")[0]
        return simple_name + "_"+ str(row["sample_id"]) +"_"+ row["primer_panel"]
    elif(execution_mode=='fastq'):
        return str(row["sample_id"]) +"_"+ row["primer_panel"] 

df["internal_pipeline_id"] = df.apply(lambda row: genPipelineInternalId(row,execution_mode), axis=1)
duplicated_bool = df["internal_pipeline_id"].duplicated(keep=False)
REQUIRED_HEADERS.append('internal_pipeline_id')
duplicated_pipe_ids = df[REQUIRED_HEADERS].loc[duplicated_bool]
try:
    assert(len(duplicated_pipe_ids)==0)
    print("    > PASS")

except(AssertionError):
    print("ERROR: The following columns will generate duplicated internal pipeline ids")
    print(duplicated_pipe_ids)
    ERRORS_FOUND +=1

# check if only valid values for primer panel were provided
print("@ checking if primer panel values provided are valid...")
def isPanelValid(row):
    if row["primer_panel"] in REQUIRED_PANELS:
        return True
    else:
        return False
isvalid_bools = df.apply(isPanelValid, axis=1)
invalid_ppnls = df.loc[~isvalid_bools]

try:
    assert(len(invalid_ppnls)==0)
    print("    > PASS")

except(AssertionError):
    print(f"ERROR: The following rows have non valid primer names as set by {panel_set_flpath}")
    print(invalid_ppnls[REQUIRED_HEADERS])
    ERRORS_FOUND +=1

try:
    assert(ERRORS_FOUND == 0)
    print(":: DONE ::")

except(AssertionError):
    print("ERRORS WERE FOUND AT IRODS MANIFEST!")
    print("FIX IT AND TRY AGAIN")
    exit(1)
