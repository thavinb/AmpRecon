#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser
from os import getcwd
from os import system

'''
--- renameSamplesCram.py ------------

This script get a sample name as output from samtools split and renames to keep
the information.
It gets sample names, index and barcode relationship from [run_id]_manifest.csv.

Assumes the names from samtools always follow the pattern
            [run_id]_*_*#[index].cram

From there it renames to
            [sample_id]#[index]_.cram

It double check the consistency between barcodes, index and sample name by comparing
the barcodes-index relationship provided at 'decode.bam.metrics'
and the barcodes-sample_name provided at [run_id]_manifest.csv
'''

# --- FUNCTIONS ---------------------------------------------------------------
def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", "-m",
                        help="path to [run_id]_manifest.csv")
    parser.add_argument("--bam_metrics", "-bm",
                        help="path to bam_metrics file")
    parser.add_argument('--cram_files', '-cfs', nargs='+', default=[])
    parser.add_argument('--run_id', default=None)
    parser.add_argument('--lane', default=None)
    return parser.parse_args()

def loadManifestDf(mnf_path):
    """
    load a dataframe containing information of a [run_id]_manifest.csv

    INPUT
    -----
    mnf_path:<str>
            path to a [run_id]_manifest.csv
    """
    # Here we will count the number of lines which should be ignored before
    # loading the manifest dataframe
    with open(mnf_path, 'r') as mnf_f:
        lines2skip = 0
        for line in mnf_f:
            if line.startswith("lims_id"):
                break
            else:
                lines2skip += 1

    return pd.read_csv(mnf_path, skiprows=lines2skip)

def parseBamMetrics(bam_metrics_path):
    """
    Load SampleSheet dataframe
    """
    def __intfy(row):
        try:
            return int(row)
        except:
            return None
    # Skip manifest header
    keys = None
    dct_lst = []
    with open(bam_metrics_path, 'r') as metrics_f:
        for line in metrics_f:
            # ignore lines
            if line.startswith("#") or line.startswith("\n"):
                continue
            # headers line
            if line.startswith("BARCODE"):
                hline=line.strip("\n").split("\t")
                keys = hline
                continue
            # data lines
            else:
                dct = {}
                dline=line.strip("\n").split("\t")
                for i, k_i in enumerate(keys):
                    dct[k_i] = dline[i]
                dct_lst.append(dct)
    metrics_df = pd.DataFrame(dct_lst)
    metrics_df["BARCODE_NAME"] = metrics_df["BARCODE_NAME"].apply(__intfy)
    return metrics_df

def getNewName(cram_fl, mnf_df, run_id, lane):
    basename = cram_fl.split('/')[-1].split('.')[0]
    cram_idx = int(basename.split('#')[-1])
    # if no idx is found, give a warning and
    try:
        sample_name = mnf_df.loc[mnf_df["index"]==cram_idx]["lims_id"].values[0]
    except(IndexError):
        print(f"WARN: no sample name for ${cram_fl}, it will be ignored")
        return None
    new_name = f"{run_id}_{lane}#{cram_idx}_{sample_name}-.cram"
    return new_name

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # parse inputs
    args = parse_args()

    # check if cram files were provided
    assert(len(args.cram_files)>0), "ERROR: no cram files were provided"

    # load bam_metrics data
    metrics_df = parseBamMetrics(args.bam_metrics)

    #metrics_df = parseBamMetrics(args.bam_metrics)

    # gets sample names, index and barcode relationship from
    # [run_id]_manifest.csv.

    mnf_df = loadManifestDf(args.manifest)
    mnf_df[["lims_id","index","barcode_sequence"]]

    # joint data on same df
    merged_df= metrics_df[["BARCODE","BARCODE_NAME"]].merge(
                mnf_df[["lims_id","index","barcode_sequence"]],
                left_on="BARCODE_NAME",
                right_on="index"
                )
    # double check if barcodes from metrics and manifest are the same
    common_df = merged_df.loc[merged_df["BARCODE"] == merged_df["barcode_sequence"]]
    errm2 =f"ERROR: barcodes don't match between {args.manifest} and {args.bam_metrics}"
    assert(len(common_df) == len(merged_df)), errm2

    # rename files
    for cram_fl in args.cram_files:
        new_name = getNewName(cram_fl, mnf_df, args.run_id, args.lane)
        # ignore cram file with no index found
        if new_name==None:
            continue
        system(f"cp --preserve=links {cram_fl} {getcwd()}/{new_name}")

