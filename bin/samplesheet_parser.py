#!/usr/bin/env python3

import logging
import re
import string
from argparse import ArgumentParser
from collections import OrderedDict
from csv import DictReader, DictWriter
from os.path import dirname, join, split
from typing import Optional
import os


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "samplesheet_path",
        help="Path to SampleSheet.csv. Preferably exists inside a run folder matching [0-9]{6}",
    )

    parser.add_argument(
        "--output_file",
        "-o",
        help="The output path to save the manifest"
    )
    return parser.parse_args()


class SampleSheetParser:
    class MissingValueError(ValueError):
        pass

    def __init__(self, samplesheet_path: str, output_file, run: Optional[int] = None):
        self.path = samplesheet_path
        self.run = run
        self.output_file = output_file
        self.reader = None
        self.writer = None
        self.well_iterator = self._well_iterator()
        self._fileio = None
        self._assay_column = None
        self._first_assay_id = None
        self._assay_iterator = 0

    def __repr__(self):
        """Shows how the class gets represented in the Pycharm debugger, for example"""
        return f"{self.__class__.__name__}: {self.path}"

    def _prepare(self):
        """Read file if it's not read yet, otherwise reset to top minus header"""
        if self._fileio is None:
            raise IOError("file not opened yet, run object inside a context manager")

        # Skip manifest header
        for row in self._fileio:
            if row.lower().startswith("[data]"):
                break

        # Get csv headers
        samplesheet_columns = next(self._fileio).lower().strip().split(",")
        manifest_columns = [
            "lims_id",
            "sims_id",
            "index",
            "assay",
            "barcode_sequence"
        ]

        self.reader = DictReader(
            self._fileio, fieldnames=samplesheet_columns, dialect="excel"
        )
        self.writer = DictWriter(self._manifest_io, fieldnames=manifest_columns)
        self.writer.writeheader()

    def _get_lims_id(self, line: OrderedDict):
        lims_id = line.get("sample_name")
        if lims_id:
            if ":" in lims_id:
                lims_id = lims_id.split(":")[1]
            elif "_" in lims_id:
                lims_id = lims_id.split("_")[0]
        else:
            raise self.MissingValueError("No sample_name present, can't assign lims id")

        return lims_id

    def _get_sims_id(self, line: OrderedDict, lims_id: str):
        sims_id = line.get("sample_id", lims_id)
        if sims_id is not None:
            if ":" in sims_id:
                sims_id = sims_id.split(":")[1]
            elif "_" in sims_id:
                sims_id = sims_id.split("_")[0]

        return sims_id

    def _get_assay(self, line: OrderedDict, lims_id: str):
        assay = None
        match = None
        assay_options = OrderedDict(
            # TODO: this needs to be taken from the panel_settings, not hardcoded.
            grc1="PFA_GRC1_v1.0", grc2="PFA_GRC2_v1.0", spec="PFA_Spec" 
        )
        grc = re.compile(r"(grc[12])|(sp(?:ec)?)", re.IGNORECASE)

        if self._first_assay_id is None:
            self._first_assay_id = lims_id
        elif self._first_assay_id == lims_id:
            self._assay_iterator += 1

        # assay can be literally everywhere apparently.. so we go search for it in every column until found
        # but once we find it we'll just reuse that column
        if self._assay_column is None:
            for header, value in line.items():
                match = re.search(grc, value)
                if match:
                    self._assay_column = header
                    break

            # apparently it's ok if assay isn't given, they always come in order of grc1, grc2, spec
            assay = assay_options[list(assay_options.keys())[self._assay_iterator]]

        else:
            match = re.search(grc, line[self._assay_column])

        if match:
            group = match.group(0).lower()
            if group.startswith("sp"):
                group = "spec"
            assay = assay_options[group]

        return assay

    def _get_barcode_sequence(self, line: OrderedDict):
        index_1 = line.get("index")
        index_2 = line.get("index2")
        if not all((index_1, index_2)):
            raise self.MissingValueError(
                "Can't find columns 'index' and 'index2' to infer barcode sequence"
            )
        return f"{index_1}-{index_2}"

    def _check_required_col(self, line:OrderedDict):
        '''
        Check if required columns are present on the samplesheet.
        '''
        # those columns are expected to be on the samplesheet
        REQUIRED_COLUMNS = ["index", "index2", "sample_name",
                            "sample_id"]
        # some of them may have an alternative name
        OPT_COLS = {"well":"sample_well",
                    "plate":"sample_plate"}

        # index and index2 are used to produce the "barcode_sequence"
        # sample_name is used to fill lims_id
        # TODO: it looks for panel names everywhere and the columns are hardcoded
        #       at _get_assay()
        missing_cols = []
        for col in REQUIRED_COLUMNS:
            if col not in line:
                # check if the optional names are available
                if col in OPT_COLS:
                    alt_nm = OPT_COLS[col]
                    if alt_nm in line:
                        continue
                missing_cols.append(col)

        if len(missing_cols) > 0:
            raise Exception(f"The following columns are missing: {missing_cols}")

    def process(self):
        """Run all validation functions"""
        for index, line in enumerate(self.reader, 1):
            # strip trailing newlines from each value
            for column in line:
                if line[column] is not None:
                    line[column] = line[column].strip()

            # check if required columns are present
            self._check_required_col(line)

            # parse values
            try:
                lims_id = self._get_lims_id(line)
                sims_id = self._get_sims_id(line, lims_id)
                assay = self._get_assay(line, lims_id)
                barcode_sequence = self._get_barcode_sequence(line)
            except self.MissingValueError as err:
                logging.error(str(err))
                exit(1)

            row = dict(
                lims_id=lims_id,
                sims_id=sims_id,
                index=index,
                assay=assay,
                barcode_sequence=barcode_sequence,
            )
            self.writer.writerow(row)

    def _well_iterator(self):
        """A generator that returns A1 - H12 as called"""
        for character in string.ascii_uppercase[:8]:  # A - H
            for number in range(1, 13):  # 1 - 12
                yield f"{character}{number}"

    def __enter__(self):
        """Makes the class function like a context manager"""
        self._fileio = open(self.path)
        self._manifest_io = open(self.output_file, "w")
        self._prepare()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the context manager and the file object"""
        self._fileio.close()
        self._manifest_io.close()


if __name__ == "__main__":
    args = parse_args()
    with SampleSheetParser(args.samplesheet_path, args.output_file, None) as s:
        s.process()
