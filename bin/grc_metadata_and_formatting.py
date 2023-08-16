#!/usr/bin/env python3

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
    "--grc1_file",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--grc2_file",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--barcodes_file",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--output_file_grc1",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--output_file_grc2",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--output_file_barcodes",
    type=str,
    required=True,
    help=".",
)

parser.add_argument(
    "--config",
    type=str,
    required=True,
    help=".",
)

args = parser.parse_args()

config = json.load(args.config)["metadata"]


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


def format_file(grc_data_structure, column_names, add_coil):
    dictionary = {}
    for id, row in grc_data_structure.items():
        new_column_names = []
        # Rename certain columns
        for column in row.keys():
            if column in column_names.keys():
                column = column_names.get(column)
            new_column_names.append(column)
        updated_row = dict(zip(new_column_names, list(row.values())))
        if add_coil == True:
            # Add COIL column with no data if GRC1
            updated_row["COIL"] = "N/A"
        dictionary[id] = updated_row
    return dictionary


def write_output_file(out_file_name, column_order, output_data):
    # Write output data to file in correct order
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

# Column names for the GRCs
grc_column_names = {
    "partner_sample_id": "External_ID",
    "collection_date": "Date_of_Collection",
    "collection_location": "Location",
    "study": "Study",
    "collection_country": "Country",
}

# ------------------------------------GRC1 PROCESSING------------------------------------

# Read in GRC1
grc1 = read_file(args.grc1_file, "ID")

# Add metadata to GRC1
metadata_grc1 = add_metadata(manifest, grc1)

# Format GRC1 - change column names and add COIL column.
formatted_grc1 = format_file(metadata_grc1, grc_column_names, True)

# Write GRC1 and associated metadata to output file in correct order
grc1_column_order = config["grc1_column_order"]

write_output_file(args.output_file_grc1, grc1_column_order, formatted_grc1)

# ------------------------------------GRC2 PROCESSING------------------------------------

# Read in GRC2
grc2 = read_file(args.grc2_file, "ID")

# Add metadata to GRC2
metadata_grc2 = add_metadata(manifest, grc2)

# Format GRC2 - change column names and add COIL column.
formatted_grc2 = format_file(metadata_grc2, grc_column_names, False)

# Write GRC2 and associated metadata to output file in correct order
grc2_column_order = config["grc2_column_order"]

write_output_file(args.output_file_grc2, grc2_column_order, formatted_grc2)

# ------------------------------------BARCODES PROCESSING------------------------------------


# Read in Barcodes files
barcodes = read_file(args.barcodes_file, "ID")

# Add metadata to Barcodes
metadata_barcodes = add_metadata(manifest, barcodes)

# Format Barcodes data - change column names
barcode_column_names = {
    "partner_sample_id": "ssID",
    "ID": "MalGEN_ID",
    "study": "Study",
}
formatted_barcodes = format_file(metadata_barcodes, barcode_column_names, False)

# Write barcodes and associated metadata to output file in correct order
barcodes_column_order = config["barcodes_column_order"]

write_output_file(args.output_file_barcodes, barcodes_column_order, formatted_barcodes)
