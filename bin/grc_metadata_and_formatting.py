#!/usr/bin/env python3

import argparse
import csv

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
grc1_column_order = [
    "ID",
    "External_ID",
    "Date_of_Collection",
    "Location",
    "Country",
    "Study",
    "Species",
    "COIL",
    "McCOIL",
    "Barcode",
    "Kelch",
    "P23:BP",
    "PfCRT",
    "PfDHFR",
    "PfDHPS",
    "PfEXO",
    "PfMDR1",
    "PGB",
]
write_output_file(args.output_file_grc1, grc1_column_order, formatted_grc1)

# ------------------------------------GRC2 PROCESSING------------------------------------

# Read in GRC2
grc2 = read_file(args.grc2_file, "ID")

# Add metadata to GRC2
metadata_grc2 = add_metadata(manifest, grc2)

# Format GRC2 - change column names and add COIL column.
formatted_grc2 = format_file(metadata_grc2, grc_column_names, False)

# Write GRC2 and associated metadata to output file in correct order
grc2_column_order = [
    "ID",
    "Study",
    "External_ID",
    "PfCRT:72",
    "PfCRT:74",
    "PfCRT:75",
    "PfCRT:76",
    "PfCRT:93",
    "PfCRT:97",
    "PfCRT:218",
    "PfCRT:220",
    "PfCRT:271",
    "PfCRT:333",
    "PfCRT:353",
    "PfCRT:371",
    "PfDHFR:16",
    "PfDHFR:51",
    "PfDHFR:59",
    "PfDHFR:108",
    "PfDHFR:164",
    "PfDHFR:306",
    "PfDHPS:436",
    "PfDHPS:437",
    "PfDHPS:540",
    "PfDHPS:581",
    "PfDHPS:613",
    "PfEXO:415",
    "PfMDR1:86",
    "PfMDR1:184",
    "PfMDR1:1034",
    "PfMDR1:1042",
    "PfMDR1:1226",
    "PfMDR1:1246",
    "PGB:127",
    "PGB:128",
    "PGB:193",
    "PGB:326",
    "PGB:356",
    "PGB:484",
]
write_output_file(args.output_file_grc2, grc2_column_order, formatted_grc2)

# ------------------------------------BARCODES PROCESSING------------------------------------


# Read in Barcodes files
barcodes = read_file(args.barcodes_file, "ID")

# Add metadata to Barcodes
metadata_barcodes = add_metadata(manifest, barcodes)

# Format Barcodes data - change column names
barcode_column_names = {
    "ID": "ssID",
    "partner_sample_id": "MalGEN_ID",
    "study": "Study",
}
formatted_barcodes = format_file(metadata_barcodes, barcode_column_names, False)

# Write barcodes and associated metadata to output file in correct order
barcodes_column_order = [
    "ssID",
    "MalGEN_ID",
    "Study",
    "Pf3D7_02_v3:376222",
    "Pf3D7_02_v3:470013",
    "Pf3D7_03_v3:656861",
    "Pf3D7_04_v3:110442",
    "Pf3D7_04_v3:881571",
    "Pf3D7_05_v3:350933",
    "Pf3D7_05_v3:369740",
    "Pf3D7_06_v3:900278",
    "Pf3D7_07_v3:1044052",
    "Pf3D7_08_v3:1314831",
    "Pf3D7_08_v3:413067",
    "Pf3D7_09_v3:900277",
    "Pf3D7_11_v3:1018899",
    "Pf3D7_11_v3:1815412",
    "Pf3D7_13_v3:1056452",
    "Pf3D7_13_v3:1466422",
    "Pf3D7_14_v3:137622",
    "Pf3D7_14_v3:2164225",
    "Pf3D7_01_v3:145515",
    "Pf3D7_03_v3:548178",
    "Pf3D7_04_v3:1102392",
    "Pf3D7_04_v3:139051",
    "Pf3D7_04_v3:286542",
    "Pf3D7_04_v3:529500",
    "Pf3D7_05_v3:796714",
    "Pf3D7_07_v3:1256331",
    "Pf3D7_07_v3:461139",
    "Pf3D7_07_v3:619957",
    "Pf3D7_08_v3:417335",
    "Pf3D7_09_v3:163977",
    "Pf3D7_10_v3:317581",
    "Pf3D7_10_v3:336274",
    "Pf3D7_11_v3:1020397",
    "Pf3D7_11_v3:1294107",
    "Pf3D7_11_v3:1935227",
    "Pf3D7_11_v3:477922",
    "Pf3D7_12_v3:1663492",
    "Pf3D7_12_v3:2171901",
    "Pf3D7_13_v3:1233218",
    "Pf3D7_13_v3:1867630",
    "Pf3D7_13_v3:2377887",
    "Pf3D7_14_v3:2355751",
    "Pf3D7_14_v3:3046108",
    "Pf3D7_02_v3:529709",
    "Pf3D7_02_v3:714480",
    "Pf3D7_03_v3:155697",
    "Pf3D7_04_v3:1037656",
    "Pf3D7_04_v3:648101",
    "Pf3D7_05_v3:1204155",
    "Pf3D7_06_v3:1282691",
    "Pf3D7_06_v3:1289212",
    "Pf3D7_07_v3:1066698",
    "Pf3D7_07_v3:1213486",
    "Pf3D7_07_v3:704373",
    "Pf3D7_08_v3:1313202",
    "Pf3D7_08_v3:339406",
    "Pf3D7_08_v3:701557",
    "Pf3D7_09_v3:452690",
    "Pf3D7_09_v3:599655",
    "Pf3D7_10_v3:1383789",
    "Pf3D7_10_v3:1385894",
    "Pf3D7_11_v3:1006911",
    "Pf3D7_11_v3:1295068",
    "Pf3D7_11_v3:1802201",
    "Pf3D7_12_v3:1667593",
    "Pf3D7_12_v3:1934745",
    "Pf3D7_12_v3:858501",
    "Pf3D7_13_v3:1419519",
    "Pf3D7_13_v3:159086",
    "Pf3D7_13_v3:2161975",
    "Pf3D7_13_v3:2573828",
    "Pf3D7_13_v3:388365",
    "Pf3D7_14_v3:2625887",
    "Pf3D7_14_v3:3126219",
    "Pf3D7_14_v3:438592",
    "Pf3D7_01_v3:179347",
    "Pf3D7_01_v3:180554",
    "Pf3D7_01_v3:283144",
    "Pf3D7_01_v3:535211",
    "Pf3D7_02_v3:839620",
    "Pf3D7_04_v3:426436",
    "Pf3D7_04_v3:531138",
    "Pf3D7_04_v3:891732",
    "Pf3D7_05_v3:172801",
    "Pf3D7_06_v3:574938",
    "Pf3D7_07_v3:1308383",
    "Pf3D7_07_v3:1358910",
    "Pf3D7_07_v3:1359218",
    "Pf3D7_07_v3:635985",
    "Pf3D7_08_v3:1056829",
    "Pf3D7_08_v3:150033",
    "Pf3D7_08_v3:399774",
    "Pf3D7_09_v3:1379145",
    "Pf3D7_10_v3:1386850",
    "Pf3D7_11_v3:1935031",
    "Pf3D7_11_v3:408668",
    "Pf3D7_11_v3:828596",
    "Pf3D7_12_v3:857245",
    "Pf3D7_14_v3:107014",
    "Pf3D7_14_v3:1757603",
    "Pf3D7_14_v3:2733656",
]
write_output_file(args.output_file_barcodes, barcodes_column_order, formatted_barcodes)
