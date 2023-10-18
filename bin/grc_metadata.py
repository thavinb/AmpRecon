#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

import argparse
import csv
import json

parser = argparse.ArgumentParser()

parser.add_argument(
    "--manifest_file",
    "-m",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--grc_file",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--output_file_grc",
    type=str,
    required=True,
    help=".",
)

args = parser.parse_args()

def read_file(input_file, key_column):
    dictionary = {}
    in_file = open(input_file, newline="")
    for row in csv.DictReader(in_file, delimiter="\t"):
        key = row.get(key_column)
        dictionary[key] = dict(row)
    in_file.close()
    return dictionary


def add_metadata(manifest, grc):
    dictionary = {}
    for id, row in grc.items():
        try:
            matching_metadata = manifest.get(id)
        except IndexError:
            raise Exception(f"Sample {id} is missing from the manifest.")
        row.update(matching_metadata)
        dictionary[id] = row
    return dictionary


def write_output_file(out_file_name, column_order, output_data):
    # Write output data to file
    output_file = open(out_file_name, "w")
    output_file_writer = csv.DictWriter(
        output_file, delimiter="\t", extrasaction="ignore", fieldnames=column_order
    )
    output_file_writer.writeheader()
    output_data_rows = output_data.values()
    output_file_writer.writerows(output_data_rows)
    output_file.close()


# Read in the manifest file
manifest = read_file(args.manifest_file, "sample_id")

# ------------------------------------GRC PROCESSING------------------------------------

# Read in GRC
grc = read_file(args.grc_file, "ID")

# Add metadata to GRC
metadata_grc = add_metadata(manifest, grc)
column_order = list(metadata_grc.values())[0].keys()

# Write GRC1 and associated metadata to output file
write_output_file(args.output_file_grc, column_order, metadata_grc)