#!/usr/bin/env python3

from ast import Assert
import pandas as pd
from sys import argv, exit

# load manifest
mnf_flpath = argv[1]
pannel_set_flpath = argv[2]
df = pd.read_csv(mnf_flpath, sep='\t')

# Set Global requirements
REQUIRED_HEADERS = ["sample_id", "primer_panel","irods_path"]
REQUIRED_FORMAT = ".cram"

pnl_df = pd.read_csv(pannel_set_flpath, sep=',')
REQUIRED_PANNELS = pnl_df["pannel_name"].to_list()
ERRORS_FOUND = 0

# check if required headers are present
print("@ checking column headers...")
columns_lst = df.columns.tolist()
found_headers = []
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

invalid_format_idx = df['irods_path'].apply(isItRequiredFormat)
non_cram_rows = df.loc[~invalid_format_idx]
try:
    assert(len(non_cram_rows) == 0)
    print("    > PASS")
except:
    print("ERROR: Non cram files were requested at rows")
    print(f"{non_cram_rows[REQUIRED_HEADERS]}")
    ERRORS_FOUND +=1

# check if WG_LANE_SAMPLE_ID is unique
print("@ check for pipeline internal ids uniqueness...")
def genPipelineInternalId(row):
    simple_name = row['irods_path'].split("/")[-1].split(".")[0]
    return simple_name + "_"+ str(row["sample_id"]) +"_"+ row["primer_panel"]

df["internal_pipeline_id"] = df.apply(genPipelineInternalId, axis=1)
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

# check if only valid values for primer pannel were provided
print("@ checking if primer panel values provided are valid...")
def isPannelValid(row):
    if row["primer_panel"] in REQUIRED_PANNELS:
        return True
    else:
        return False
isvalid_bools = df.apply(isPannelValid, axis=1)
invalid_ppnls = df.loc[~isvalid_bools]

try:
    assert(len(invalid_ppnls)==0)
    print("    > PASS")

except(AssertionError):
    print(f"ERROR: The following rows have non valid primer names as set by {pannel_set_flpath}")
    print(invalid_ppnls[REQUIRED_HEADERS])
    ERRORS_FOUND +=1

try:
    assert(ERRORS_FOUND == 0)
    print(":: DONE ::")

except(AssertionError):
    print("ERRORS WERE FOUND!")
    exit(1)