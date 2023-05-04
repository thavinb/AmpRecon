#!/usr/bin/env python3

import argparse
import json
from glob import glob
import os
from multiprocessing import Pool
import csv
from collections import defaultdict


class NoConfigError(Exception):
    def __init__(
        self,
        message="""
        You must specify a path to a valid config file.
    """,
    ):
        self.message = message
        super().__init__(self.message)


class NoGenotypeFilesError(Exception):
    def __init__(
        self,
        message="""
        You must specify input genotype files.
    """,
    ):
        self.message = message
        super().__init__(self.message)


class NoBarcodeRefError(Exception):
    def __init__(
        self,
        message="""
        Please provide a config file containing a barcode_ref.
    """,
    ):
        self.message = message
        super().__init__(self.message)


class InvalidNCPUsRequested(Exception):
    def __init__(
        self,
        message="""
        You cannot request less than 1 CPU.
    """,
    ):
        self.message = message
        super().__init__(self.message)


class Barcode:
    def __init__(self, d_genotypes: dict, barcode_ref: dict) -> None:
        self.barcode_ref = barcode_ref
        self.barcode_indicies = sorted(list(self.barcode_ref.keys()))
        self.d_genotypes = d_genotypes
        self.barcode = self._build_barcode()

    def _build_barcode(self) -> str:
        """
        1) iterate through sorted barcode indices
        2) at each position get barcode ref (chrom and locus)
        3) get genotype for position from genotype file
        4) Check if genotype missing or het and add to output string accordingly
        5) return barcode string
        """
        out_barcode = ""

        for index in self.barcode_indicies:
            d_barcode_pos = self.barcode_ref[index]
            gen = self.d_genotypes.get(d_barcode_pos["Chromosome"], {}).get(
                str(d_barcode_pos["Locus"])
            )
            if not gen or gen == "-":
                gen = "X"
            elif len(gen.split(",")) > 1:
                gen = "N"
            out_barcode += gen
        return out_barcode


def output_df(records, output_file: str, sep="\t") -> None:
    """
    Write out barcodes to tsv file with a location provided by the user.
    """
    with open(output_file, "w+") as outBarcodes:
        writer = None
        for row in records:
            if not writer:
                fields = list(row.keys())
                writer = csv.DictWriter(outBarcodes, fieldnames=fields, delimiter=sep)
                writer.writeheader()
            writer.writerow(row)
        outBarcodes.close()

def output_barcode_split_out(barcode_ref_dct, records, output_flnm:str, sep="\t"):
    """
    Write a second barcode file in which there is on column per barcode base and columns
    are named as {Chromosome}:{Locus} 
    """
    # get col names and positions
    cols_dct = {"cols_bar":{}}
    for pos_i in barcode_ref_dct.keys():
        chrm_i = barcode_ref_dct[pos_i]["Chromosome"]
        locs_i = barcode_ref_dct[pos_i]["Locus"]
        cols_dct["cols_bar"][int(pos_i)]= f"{chrm_i}:{locs_i}"

    # get numbers of pos expected
    n_pos = len(cols_dct["cols_bar"].keys())

    # --- sanity check ---
    # assert we have a continuous list of integers
    assert(min(cols_dct["cols_bar"].keys()) == 1) # it shoudl start with one
    max_pos = max(cols_dct["cols_bar"].keys())
    assert(set(cols_dct["cols_bar"].keys()) == set(range(1, max_pos+1))), "Not a continuous list of integers for the barcode positions"

    # get samples data
    out_sample_dct = []
    for sample_j in records:
        # record id
        smpl_dct = {"ID":sample_j["ID"]}
        
        # get bases from barcode and associate 
        # with the respective column
        barcode = sample_j["Barcode"]
        assert(len(barcode) == n_pos)
    
        for i, base_i in enumerate(barcode):
            colnm = cols_dct["cols_bar"][i+1]
            smpl_dct[colnm] = base_i
    
        out_sample_dct.append(smpl_dct)

    # write tsv file
    with open(output_flnm,"w") as out_file:
        writer = None
        for row in out_sample_dct:
            if not writer:
                fields = list(row.keys())
                writer = csv.DictWriter(out_file, fieldnames=fields, delimiter=sep)
                writer.writeheader()
            writer.writerow(row)
        out_file.close()

def read_genotype_file(genotype_file_path: str, deconstruct_barcode_ref: dict) -> dict:
    """
    Read in genotype file, extract records which are in the barcode ref into dictionary.

    Output dictionary structured as: {CHROM:{LOC:GEN}}
    """
    with open(genotype_file_path) as genotypefile:
        reader = csv.DictReader(genotypefile, delimiter="\t")

        d_out = defaultdict(lambda: defaultdict(dict))

        for row in reader:
            if (
                row["Loc"] in deconstruct_barcode_ref
                and row["Chr"] in deconstruct_barcode_ref[row["Loc"]]
            ):
                d_out[row["ID"]][row["Chr"]][row["Loc"]] = row["Gen"]
        genotypefile.close()

    return dict(d_out)


def main(
    genotype_file_path: str, barcode_ref: dict, deconstruct_barcode_ref: dict
) -> list:
    """
    1) Load genotype file into dictionary with read_genotype_file
    2) Create Barcode
    3) Output dictionary with sample name and barcode for each sample in genotype file dict
    """
    out_rows = []
    d_genotype_file = read_genotype_file(genotype_file_path, deconstruct_barcode_ref)

    for sample, d_sample_genotype in d_genotype_file.items():
        barcoding = Barcode(d_sample_genotype, barcode_ref)
        out_rows.append({"ID": sample, "Barcode": barcoding.barcode})
    return out_rows


def map_main(args):
    """
    Return run of the main function with arguments
    (Implementation of multiprocessing starmap)
    """
    return main(*args)


if __name__ == "__main__":
    # initialise argument parser
    parser = argparse.ArgumentParser(
        prog="grc_barcoding.py",
        description="A package to perform barcode production based on the production amplicon pipeline",
    )

    # create argparser arguments
    parser.add_argument(
        "--genotype_files",
        help="""
        Path to input genotype file(s)
        
    """,
        nargs="+",
    )
    parser.add_argument(
        "--config",
        help="""
        Path to config json file
    """,
    )
    parser.add_argument(
        "--output_file",
        help="""
        Path to directory to output results (default: barcode_results.txt)
    """,
    )
    parser.add_argument(
        "--pbar",
        action="store_true",
        help="""
        Show a progress bar while running 
    """,
    )
    parser.add_argument(
        "--ncpus",
        type=int,
        help="""
        No. cpus to use in processing
    """,
    )

    parser.add_argument(
        "--output_file_split_out",
        help="""
        Path to directory to output results (example: barcoding_results.split_out.txt)
    """,
    )

    # get all arguments with a value and place into dictionary
    args = {k: v for k, v in vars(parser.parse_args()).items() if v}

    # read in config json file provided at command line, get speciation portion
    # if not provided throw error
    try:
        config = json.load(open(args.pop("config")))["grc_barcoding"]
    except KeyError:
        raise NoConfigError

    # get the barcode ref from the config, if not present throw an error.
    try:
        barcode_ref = config.pop("barcode_ref")
    except:
        raise NoBarcodeRefError

    # populate deconstruct barcode for keying in read_genotype_file
    deconstruct_barcode_ref = defaultdict(list)
    for d in barcode_ref.values():
        deconstruct_barcode_ref[str(d["Locus"])].append(d["Chromosome"])

    # if genotype files not provided throw error
    if not args["genotype_files"]:
        raise NoGenotypeFilesError

    # get list of genotype file paths and expand any glob string file paths
    genotype_files = []
    for path in args["genotype_files"]:
        genotype_files.extend(glob(path))

    # get optional arguments if provided, if not return default value
    output_file = os.path.abspath(args.get("output_file", "./barcoding_results.txt"))
    ncpus = args.get("ncpus", 1)
    pbar = args.get("pbar", False)
    output_file_2 = os.path.abspath(args.get("output_file_split_out", "./barcoding_results.split_out.txt"))
    # add in logic to condionally attempt to load tqdm for pbar
    # if tqdm not present continue without progress bar
    if pbar:
        try:
            from tqdm import tqdm
        except ModuleNotFoundError:
            print(
                """
            TQDM not installed, continuing with no progress bar.
            To use progress bar please install TQDM: https://github.com/tqdm/tqdm
            """
            )
            pbar = False

    # check if ncpus is one, if so run in serial (bypass multiprocessing)
    if ncpus == 1:
        # check if user has asked for pbar, if  so create iterator with tqdm
        # otherwise return iterator as is
        if pbar:
            iterator = tqdm(genotype_files)
        else:
            iterator = genotype_files

        # initialise all_samples list to collect data for each genotype file
        all_samples_out = []
        for genotype_file in iterator:
            # run a single genotype file through main
            out = main(genotype_file, barcode_ref, deconstruct_barcode_ref)
            # append output to all_samples_out for output
            all_samples_out.append(out)
    elif ncpus > 1:
        # if ncpus >1 set up argument tuples per sample
        all_args = [
            (genotype_file, barcode_ref, deconstruct_barcode_ref)
            for genotype_file in genotype_files
        ]

        # if pbar requested create iterator with tqdm, otherwise don't
        with Pool(ncpus) as p:
            if pbar:
                all_samples_out = list(
                    tqdm(p.imap(map_main, all_args), total=len(genotype_files))
                )
            else:
                all_samples_out = list(p.imap(map_main, all_args))
    else:
        raise InvalidNCPUsRequested

    # flatten nested list from output
    all_samples_out = [i for sublist in all_samples_out for i in sublist]
    # pass flattened list to output_df
    output_df(all_samples_out, output_file)

    output_barcode_split_out(barcode_ref, all_samples_out,output_file_2)