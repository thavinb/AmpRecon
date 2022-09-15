#!/usr/bin/env python3

import logging
import os
from argparse import ArgumentParser
from collections import OrderedDict
from csv import DictReader


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", "-m", help="path to manifest to validate")
    return parser.parse_args()


class AmpliconManifestValidator:
    class InvalidColumnsError(ValueError):
        pass

    class InvalidValueError(ValueError):
        pass

    class EmptyValuesError(ValueError):
        pass

    def __init__(self, manifest_path: str):
        self.path = manifest_path
        self.base = os.path.basename(self.path)
        self.reader = None
        self.assay_set = set()
        self.valid_assay_values = ["PFA_GRC1_v1.0", "PFA_GRC2_v1.0", "PFA_Spec"]
        self._fileio = None

    def __repr__(self):
        """Shows how the class gets represented in the Pycharm debugger, for example"""
        return f"{self.__class__.__name__}: {self.path}"

    def validate(self):
        """Run all validation functions"""
        self._prepare()
        try:
            self._validate_integrity()
            for line in self.reader:
                for key in line:
                    if line[key] is not None:
                        line[key] = line[key].strip()
                self._check_empty_or_na(line)
                self._validate_assay_column(line)
                self._validate_index_column(line)
                self._validate_barcode(line)

            if len(self.assay_set) < 3:
                raise self.InvalidValueError(
                    f"assay column does not contain all assayerences: {', '.join(self.valid_assay_values)}"
                )
        except (
            self.InvalidValueError,
            self.InvalidColumnsError,
            self.EmptyValuesError,
        ) as err:
            logging.error(str(err))
            exit(1)

    def __enter__(self):
        """Makes the class function like a context manager"""
        self._fileio = open(self.path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the context manager and the file object"""
        self._fileio.close()

    def _prepare(self):
        """Read file if it's not read yet, otherwise reset to top minus header"""
        if self._fileio is None:
            raise IOError("file not opened yet, run object inside a context manager")

        # Skip manifest header
        for row in self._fileio:
            print(row)
            if row.lower().startswith("[data]"):
                break

        # get manifest columns
        #print(next(self._fileio))#.lower().strip().split(","))
        columns = next(self._fileio).lower().strip().split(",")

        self.reader = DictReader(self._fileio, fieldnames=columns)

    def _validate_integrity(self):
        """Validates integrity of the manifest, like columns, line length, and empty/NA values"""
        valid_columns = [
            "lims_id",
            "sims_id",
            "index",
            "assay",
            "barcode_sequence",
            "well",
            "plate",
        ]

        actual_columns = self.reader.fieldnames
        if actual_columns != valid_columns:
            raise self.InvalidColumnsError(
                f"{self.base} - unexpected column names. "
                f"Expected {', '.join(valid_columns)} - got {', '.join(actual_columns)}"
            )

    def _check_empty_or_na(self, line: OrderedDict):
        for key, value in line.items():
            if not value:
                raise self.EmptyValuesError(
                    f"{self.base} has empty values: {key, value} - {line}"
                )
            elif value.upper() == "NA":
                raise self.InvalidValueError(
                    f"{self.base} has NA values: {key, value} - {line}"
                )

    def _validate_assay_column(self, line: OrderedDict):
        """Checks whether the value in the assay column matches expected values"""
        try:
            test = line["assay"]
        except KeyError:
            # Something is wrong with the file, should be picked up by other method
            # This simply assumes that validate_integrity hasn't been called yet
            self._validate_integrity()
        else:
            if test not in self.valid_assay_values:
                raise self.InvalidValueError(
                    f"{self.base} - Invalid value in assay column. "
                    f"Expected one of: {', '.join(self.valid_assay_values)} - got {test}"
                )
            else:
                self.assay_set.add(test)

    def _validate_index_column(self, line: OrderedDict):
        """Checks whether the index is a number"""
        try:
            int(line["index"])
        except KeyError:
            # Something is wrong with the file, should be picked up by other method
            # This simply assumes that validate_integrity hasn't been called yet
            self._validate_integrity()
        except ValueError:
            # index is not a number, int threw an error
            raise self.InvalidValueError(
                f"{self.base} - Invalid value in index column. Expected integer - got {line['index']}"
            )

    def _validate_barcode(self, line: OrderedDict):
        """Checks that the barcode sequence follows two DNA sequences separated by a hyphen format"""
        nucleotides = {"A", "C", "G", "T"}
        try:
            barcode = line["barcode_sequence"]
        except KeyError:
            # Something is wrong with the file, should be picked up by other method
            # This simply assumes that validate_integrity hasn't been called yet
            self._validate_integrity()
        else:
            split = barcode.split("-")

            # barcode has no hyphen separator
            if split[0] == barcode:
                raise self.InvalidValueError(
                    f"{self.base} - Invalid value in barcode sequence column - missing separator. {barcode}"
                )

            for seq in split:  # We should have 2 short sequences now
                # Check if all items from seq are in nucleotides
                if not set(seq).issubset(nucleotides):
                    raise self.InvalidValueError(
                        f"{self.base} - Invalid value in barcode sequence column. {barcode}"
                    )


if __name__ == "__main__":
    args = parse_args()
    with AmpliconManifestValidator(args.manifest) as m:
        m.validate()

