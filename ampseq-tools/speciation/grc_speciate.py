#!/usr/bin/env python3

from speciation.speciation import Speciate
import os
import json
import argparse
from glob import glob
from pathlib import Path
from multiprocessing import Pool
import csv
from collections import defaultdict
import re


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


class NoBarcodesError(Exception):
    def __init__(
        self,
        message="""
        You must specify a path to a valid barcodes file.
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


def output_df(records, output_file: str, sep="\t") -> None:
    """
    Write out barcodes to tsv file with a location provided by the user.
    """
    with open(output_file, "w+") as outBarcodes:
        writer = None
        for row in records:
            if not writer:
                fields = list(row.keys())
                writer = csv.DictWriter(outBarcodes, fieldnames=fields, delimiter="\t")
                writer.writeheader()
            writer.writerow(row)
        outBarcodes.close()


def read_barcodes(tsv_path):
    """
    Read in barcodes file into a dictionary with key,value pairs of sampleID:Barcode
    """
    out = {}
    with open(tsv_path) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        for row in reader:
            out[row["ID"]] = row["Barcode"]
    return out


def read_genotype_file(
    genotype_file_path: str, species_ref: list, chrom_regex: str
) -> dict:
    """
    Read in a genotype file as a TSV and output a dictionary in the format of:

    {
        sampleID:{
            Amplicon:{
                Locus:{
                    Depth,
                    Genotype
                }
            }
        }
    }
    """
    with open(genotype_file_path) as genotypefile:
        reader = csv.DictReader(genotypefile, delimiter="\t")

        d_out = defaultdict(lambda: defaultdict(dict))

        samples = defaultdict(bool)

        for row in reader:
            samples[row["ID"]] = True
            if (
                re.match(chrom_regex, row["Amplicon"])
                and str(row["Loc"]) in species_ref
                and row["Depth"] != "-"
            ):
                d_out[row["ID"]][row["Amplicon"]][row["Loc"]] = {
                    "Depth": row["Depth"],
                    "Gen": row["Gen"],
                }
                samples[row["ID"]] = True
        for sample, switch in samples.items():
            if not switch:
                d_out[sample] = {}

    return dict(d_out)


def main(
    genotype_file_path: str,
    d_barcode: str,
    maf_out: str,
    config: dict,
    chrom_regex: str,
) -> dict:
    """
    1) Load genotype file from path into a dictionary
    2) Create Speciate object with genotype file, barcode and config
    3) Output tsv for MAF/Depth counts
    4) dictionary with rows for species call and atched loci outputs
    """
    d_genotype_file = read_genotype_file(
        genotype_file_path, config["species_ref"], chrom_regex
    )

    out_species_labels = []
    out_matched_loci = []

    for sample, d_sample_genotype in d_genotype_file.items():
        speciate = Speciate(d_sample_genotype, d_barcode[sample], **config)

        out_species_labels.append({"ID": sample, "Species": speciate.species_label})
        speciate.matched_loci["ID"] = sample
        out_matched_loci.append(speciate.matched_loci)

        if maf_out:
            maf_rows = [{"Position": pos} | d for pos, d in speciate.write_out.items()]
            output_df(maf_rows, maf_out, f"{sample}_maf.tsv")

    return {"species_out": out_species_labels, "matched_loci": out_matched_loci}


def map_main(args):
    """
    Return run of the main function with arguments
    (Implementation of multiprocessing starmap)
    """
    return main(*args)


if __name__ == "__main__":
    # initialise argument parser
    parser = argparse.ArgumentParser(
        prog="Amplicon Pipeline Speciation",
        description="A package to perform speciation based on the production amplicon pipeline",
    )

    # create argparser arguments
    parser.add_argument(
        "--genotype_files",
        help="""
        Input list of all genotype files (TSV format) to be run through the speciation program
        """,
    )
    parser.add_argument(
        "--barcodes_file",
        help="""
        Path to barcodes output file for querying
    """,
    )
    parser.add_argument(
        "--config",
        help="""
        Path to config json
    """,
    )
    parser.add_argument(
        "--output_file",
        help="""
        Path to output file
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
        "--output_debug_path",
        help="""
        Ouput directory to save debug files. If not provided files will not be output.
    """,
    )

    # get all arguments with a value and place into dictionary
    args = {k: v for k, v in vars(parser.parse_args()).items() if v}

    # get list of genotype files from glob string, if genotype files not provided
    # throw error
    if args["genotype_files"]:
        genotype_files = glob(args["genotype_files"])
    else:
        raise NoGenotypeFilesError

    # check for barcodes file, if present read into dataframe and transform to dict
    # for querying. If not present throw error
    if args["barcodes_file"]:
        d_barcode = read_barcodes(args["barcodes_file"])
    else:
        raise NoBarcodesError

    # read in config json file provided at command line, get speciation portion
    # if not provided throw error
    if args["config"]:
        config = json.load(open(args["config"]))["grc_speciation"]
    else:
        raise NoConfigError

    # pop the chrom regex from the config, if not present get default
    chrom_regex = config.pop("chrom_regex", "^Spec_[12]_(falciparum|vivax)$")

    # get optional arguments if provided, if not return default value
    output_file = args.get("output_file", "./species.tsv")
    debug_outdir = args.get("output_debug_path", "")
    ncpus = args.get("ncpus", 1)
    pbar = args.get("pbar", False)

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
            out = main(genotype_file, d_barcode, debug_outdir, config, chrom_regex)
            # append output to all_samples_out for output
            all_samples_out.append(out)
    elif ncpus > 1:
        # if ncpus >1 set up argument tuples per sample
        all_args = [
            (genotype_file, d_barcode, debug_outdir, config, chrom_regex)
            for genotype_file in genotype_files
        ]

        # create multiprocessing pool
        with Pool(ncpus) as p:
            # if pbar requested create iterator with tqdm, otherwise don't
            if pbar:
                all_samples_out = list(
                    tqdm(p.imap(map_main, all_args), total=len(genotype_files))
                )
            else:
                all_samples_out = list(p.imap(map_main, all_args))
    else:
        raise InvalidNCPUsRequested

    # set up lists to record rows for output dataframe
    species_out_records = []
    sample_stats_records = []

    # iterate over all genotype file outputs, extend relevant lists
    for record in all_samples_out:
        species_out_records.extend(record["species_out"])
        sample_stats_records.extend(record["matched_loci"])

    # output species info to desired destination
    output_df(species_out_records, output_file)

    # output matched loci information to outdir if requested
    if debug_outdir:
        output_df(sample_stats_records, os.path.join(debug_outdir, "matched_loci.tsv"))
