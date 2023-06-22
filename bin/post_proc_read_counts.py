#!/usr/bin/env python3
import pandas as pd
import argparse
import os.path

'''
This script 1) rename the "Rpt" column values (file_ids) to be only the lanelet, \
2) remove the "Ref_region" column and 3) reorder values by index.

PS: this is a requirement for Sanger specific analysis, by adding this we tied \
    the current file_id structure ({lanelet}_{sample_id}_{panel}) to the inner \
    works of the pipeline. If this change, this post processing need to be reviewed.
'''
# --- FUNCTIONS --- #

# rename Rpt to lanelets
def load_lanelets(row):
    return "_".join(row["Rpt"].split("_")[:2])

# get indexes
def load_indexes(row):
    return int(row["lanelet"].replace("_T","").split("_")[:2][-1].split("#")[-1])

# ------------------ #

# --- get inputs --- #
parser = argparse.ArgumentParser(
    prog="post_proc_read_couts.py",
    description="""
        This script 1) rename the "Rpt" column values (file_ids) to be only the lanelet, \
        2) remove the "Ref_region" column and 3) reorder values by index.
        """,
)

parser.add_argument(
    "-csv_in",
    help="""
    Path to read counts csv files.
""",
    required=True
)

parser.add_argument(
    "-csv_out",
    help="""
    name of the grc output (default = final.grc)
""",
    required=True
)
args = {k: v for k, v in vars(parser.parse_args()).items()}

# --------------------- #

# store input values
csv_in = args["csv_in"]
csv_out = args["csv_out"]

# --- SANITY CHECK --- #
# - file must exist
try:
    assert(os.path.isfile(csv_in) == True)
except(AssertionError):
    print(f"ERROR: {csv_in} does not exists")
    exit(1)
# -------------------

print("@ loading dataframe...")
# load dataframe
read_cnts_df = pd.read_csv(csv_in)

print("@ post processing...")
# load lanelet columns
read_cnts_df["lanelet"] = read_cnts_df.apply(load_lanelets, axis=1)

# load index columns
read_cnts_df["index"] = read_cnts_df.apply(load_indexes, axis=1)

# sort df by index, drop undesired columns and rename lanelet as Rpt
final_df = read_cnts_df.sort_values(by="index").drop(labels=["Ref_region", "Rpt","index"], axis=1)
final_df.rename(columns={"lanelet":"Rpt"}, errors="raise",inplace=True)

# assert is not empty
assert(len(final_df)>0)

# write file
print(f"@ writing {csv_out}")
my_cols_lst = list(final_df.columns)
# set order of columns ("Rpt" needs to be the first, not the last)
my_cols_lst.insert(0, my_cols_lst.pop(-1))

final_df[my_cols_lst].rename(columns={"Unnamed: 8":"", "Unnamed: 15":""}).to_csv(csv_out, index=False)
print(":: DONE ::")
