#!/usr/bin/env python

import argparse

import pandas as pd


def manifest2taglist(
    manifest: str, library_name: str, sample_name: str, study_id: str
) -> None:
    """function which generates taglist file from manifest"""
    # --- Antonio Quick Fix Alert ---------------------------------------------
    # The header is kept from the SampleSheet for reasons
    # and sometimes they vary in the number of lines. Here we will count the
    # the number of lines which should be ignored before loading the manifest
    with open(manifest, 'r') as mnf_f:
        lines2skip = 0
        for line in mnf_f:
            if line.startswith("lims_id"):
                break
            else:
                lines2skip += 1
    # -------------------------------------------------------------------------
    # read in manifest
    manifest = pd.read_csv(manifest, skiprows=lines2skip)
    # lowercase the column names
    manifest.columns = [x.lower() for x in manifest.columns]

    # rename columns for tag file
    manifest.rename(
        columns={"index": "barcode_name"},
        inplace=True,
    )

    # remove other columns
    print(manifest.columns)
    manifest = manifest[["barcode_sequence", "barcode_name"]]

    # add other columns to tag_file
    manifest["library_name"] = library_name
    manifest["sample_name"] = sample_name
    manifest["study_id"] = study_id

    # create tag_list
    tag_list = manifest
    tag_list.to_csv("tag_file.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", "-m", help="Manifest file", required=True)
    parser.add_argument("--library", "-l", help="Library name", required=True)
    parser.add_argument("--sample", "-sa", help="Sample ID", required=True)
    parser.add_argument("--study", "-st", help="Study ID", required=True)
    args = parser.parse_args()

    manifest2taglist(args.manifest, args.library, args.sample, args.study)
