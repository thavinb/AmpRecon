#!/usr/bin/env python3

#from speciation.speciation import Speciate
import os
import json
import argparse
from glob import glob
from pathlib import Path
from multiprocessing import Pool
import csv
from collections import defaultdict
import re

from collections import defaultdict

d_species_convert = {
    "falciparum":"Pf",
    "vivax":"Pv"
}

class Speciate:
    def __init__(
        self,
        d_genotype_file:dict,
        barcode:str,
        species_ref:dict,
        min_maf=0.01,
        match_threshold=0.95,
        default_species="Pf",
        default_species_barcode_coverage=51,
        d_species_convert=d_species_convert,
    ) -> None:
        #initialise variables
        self.min_maf = min_maf
        self.match_threshold = match_threshold
        self.species_ref = species_ref
        self.d_species_convert = d_species_convert
        self.write_out = {}

        #run main logic flow
        self.alleles_depth_dict = self._read_genotype_file(d_genotype_file)
        self._manipulate_maf()
        self.merged_alleles = self._merge_alleles()
        self.matched_loci = self._match_loci()
        self.species_label = self._get_species_tag(
            barcode, 
            default_species, 
            default_species_barcode_coverage
            )

    @staticmethod
    def _create_summed_depth(depth: str) -> list:
        depth = sum([int(i) for i in depth.split(",")])
        return depth

    def _read_genotype_file(
        self, 
        d_genotype_file: dict
        ) -> dict:
        """
        1) get genotype file dataframe
        2) remove rows which aren't in keep_choms list
        3) Remove rows where depth is -
        4) Create summed depth from depth column
        5) Group by amplicon, get speices, iterate through records recording alleles and depth in out
        6) for each position in out if one species not present add dict with empty entries
        """
        out = defaultdict(dict)
        
        for chrom, d_chrom in d_genotype_file.items():
            species = self.d_species_convert[chrom.split("_")[-1]]
            for loc, d_loc in d_chrom.items():
                summedDepth = self._create_summed_depth(d_loc["Depth"])
                out[loc][species] = {
                    "Allele":[i for i in list(d_loc["Gen"]) if i!=","], #flexible to column being comma separated or not
                    "DP":summedDepth
                }

        #Iterate through and apply empty rows to any position missing a species call
        out_copy = out.copy()
        for pos, d_gt in out_copy.items():
            for species in self.d_species_convert.values():
                if species not in d_gt:
                    out[pos][species] = {"Allele":[],"DP":0}
        out_copy = None
        return out

    @staticmethod
    def _calculate_maf(speciesDepth: int, totDepth: int) -> float:
        """
        Calculate Pf/Pv MAF
        """
        return 1-int(speciesDepth)/int(totDepth)

    def _manipulate_maf(self) -> None:
        """
        Iterate through each position in alleles_depth_dict
        1) get total depth for each position
        2) Check totDepth > 0
        3) if Pf DP > Pv DP use Pf DP for MAF calculation, vice versa
        4) Check if there is a species (won't be if no depth/depths equal)
        5) calculate maf
        6) If MAF < min_maf (0.01 in production), remove allele calls from other species
        """
        #will be changing alleles_depth_dict during iteration,
        #therefore iterate a copy
        iterable = self.alleles_depth_dict.copy()
        
        for pos, d_gt in iterable.items():
            species=None
            totDepth=d_gt["Pf"]["DP"]+d_gt["Pv"]["DP"]
            #only apply if there is a total depth greater than 0
            #otherwise no modification will occur to that row
            maf=0
            DP=0
            if totDepth>0:
                if d_gt["Pf"]["DP"]>d_gt["Pv"]["DP"]:
                    species="Pf"
                    DP=d_gt["Pf"]["DP"]
                elif d_gt["Pv"]["DP"]>d_gt["Pf"]["DP"]:
                    species="Pv"
                    DP=d_gt["Pv"]["DP"]
                
                if DP:
                    #calculate MAF based on whichever species depth consists of the major allele
                    maf = self._calculate_maf(DP, totDepth)
                    #if maf<threshold then remove minor species allele (whether that is Pf or Pv)
                    if maf<self.min_maf:
                        if species=="Pf":
                            self.alleles_depth_dict[pos]["Pv"]["Allele"] = [] 
                        elif species=="Pv":
                            self.alleles_depth_dict[pos]["Pf"]["Allele"] = []
                #if no species set that indicates Pf and Pv have equal depths, set maf as 0.5
                if not species:
                    maf=0.5
            
            #for logging: write Pf/Pv depth and MAF to write_out dictionary
            if pos in self.species_ref:
                self.write_out[pos] = {
                    "PfDepth":d_gt["Pf"]["DP"],
                    "PvDepth":d_gt["Pv"]["DP"],
                    "MAF":maf
                }
        #remove iterable for memory purposes
        iterable=None

    def _merge_alleles(self) -> dict:
        """
        For each position in alleles depth dict iterate through each species
        1) add all unique alleles for that position from both species alignments
           to out dict
        2) return dict with each species list sorted
        """
        out = defaultdict(set)
        for pos, d_depth in self.alleles_depth_dict.items():
            for species in self.d_species_convert.values():
                #combine unique alleles for both Pf and Pv alignments
                out[pos].update(d_depth[species]["Allele"])
                if pos in self.species_ref:
                    #log unique alleles to write_out
                    self.write_out[pos]["Call"] = ",".join(list(out[pos]))
        return {k:sorted(list(v)) for k,v in out.items()}

    def _match_loci(self) -> dict:
        """
        Calculate frequency of each species specific position called

        For each position in species ref
        1) get alleles from sample (if missing get empty list)
        2) For each species in that position from species ref add 
           to count if allele exists in sample alleles
        3) Create proportions from no. discriminating positions with a call
        4) return dict
        """
        out = defaultdict(int)
        tot=0
        #iterate through all species_ref positions
        for pos, d_alleles in self.species_ref.items():
            #get the position for merged alleles, return empty list if not present
            sample_alleles = self.merged_alleles.get(pos,[])
            if sample_alleles:
                #denominator is the total non-missing species
                #discriminating positions
                tot+=1
            #iterate through species and alleles at position
            for species, allele in d_alleles.items():
                #add to species count if allele present
                if allele in sample_alleles:
                    out[species]+=1
        #divide each value by denominator "tot"
        out = {k:v/tot for k,v in out.items()}
        return out
    
    def _get_species_tag(self, barcode:str, default_species:str, min_coverage: int) -> str:
        """
        Join together species tags which have a frequency greater than match_threshold (0.95 in production).

        If barcode has sufficient coverage (51 non missing positions in production) append 
        default species (Pf in production)
        """
        #obtain list of all species with matched loci greater than threshold
        all_species_out = [k for k,v in self.matched_loci.items() if v>=self.match_threshold]

        #calculate barcode coverage as no. non missing positions in barcode
        barcode_coverage = len([i for i in list(barcode) if i!="X"])

        #chack barcode coverage above minimum, that there are species called and that the defalt species isn't
        #already in all_species_out
        if barcode_coverage >= min_coverage and all_species_out and default_species not in all_species_out:
            #if not insert the default species to the start of the list
            all_species_out.insert(0, default_species)

        #sort all species and join by backslash
        label = "/".join(sorted(all_species_out))

        if not label:
            #if no species called, output dash
            label = "-"

        return label

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
        prog="grc_speciate.py",
        description="A package to perform speciation based on the production amplicon pipeline",
    )

    # create argparser arguments
    parser.add_argument(
        "--genotype_files",
        help="""
        Input list of all genotype files (TSV format) to be run through the speciation program
        """,
        nargs="+",
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

    # if genotype files not provided throw error
    if not args["genotype_files"]:
        raise NoGenotypeFilesError

    # get list of genotype file paths and expand any glob string file paths
    genotype_files = []
    for path in args["genotype_files"]:
        genotype_files.extend(glob(path))

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
