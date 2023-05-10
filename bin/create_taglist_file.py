#!/usr/bin/env python3

import argparse
import pandas as pd


def __get_library_name(row):
    """
    assemble library name as [sample_id]_[primer_panel]
    """
    return f"{row['sample_id']}_{row['primer_panel']}"


def manifest2taglist(manifest: str, study_id: str) -> None:
    """function which generates taglist file from manifest"""
    # read in manifest
    manifest_df = pd.read_csv(manifest)  # , skiprows=lines2skip)
    # lowercase the column names
    manifest_df.columns = [x.lower() for x in manifest_df.columns]

    # rename columns for tag file
    manifest_df.rename(
        columns={"barcode_number": "barcode_name"},
        inplace=True,
    )

    # remove other columns
    tag_list_df = manifest_df[["barcode_sequence", "barcode_name"]]

    # add other columns to tag_file
    tag_list_df["library_name"] = manifest_df.apply(__get_library_name, axis=1)
    tag_list_df["sample_name"] = manifest_df["sample_id"]
    tag_list_df["study_id"] = study_id

    # create tag_list
    tag_list_df.to_csv("tag_file.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", "-m", help="Manifest file", required=True)
    parser.add_argument("--study", "-st", help="Study ID", required=True)
    args = parser.parse_args()

    manifest2taglist(args.manifest, args.study)
