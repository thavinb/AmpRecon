import pandas as pd
from os import getcwd
from os import system
from argparse import ArgumentParser


'''
--- getSampleRef.py ------------

This script get a given sample, check on the [run_id]_manifest.csv which pannel
should be used (currently, grc1, grc2 and spec) it should be aligned, and create
symlinks of the respective files to the current working directory.

PS: This script is expected to be used on the nextflow context.
'''

# --- FUNCTIONS ---------------------------------------------------------------

def loadManifestDf(mnf_path):
    '''
    load a dataframe containing information of a [run_id]_manifest.csv

    INPUT
    -----
    mnf_path:<str>
            path to a [run_id]_manifest.csv
    '''
    # Here we will count the number of lines which should be ignored before
    # loading the manifest dataframe
    with open(mnf_path, 'r') as mnf_f:
        lines2skip = 0
        for line in mnf_f:
            if line.startswith("lims_id"):
                break
            else:
                lines2skip += 1

    return pd.read_csv(mnf_path, skiprows=lines2skip)

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--manifest", "-m",
                        help="path to [run_id]_manifest.csv")
    parser.add_argument("--cram_fl", "-cfl",
                        help="path to original cram file")
    parser.add_argument('--ref_database', '-rdb',
                        help="")
    return parser.parse_args()

# -----------------------------------------------------------------------------
# --- MAIN --------------------------------------------------------------------

args = parse_args()

manifest = args.manifest
cram_fl = args.cram_fl
ref_database = args.ref_database

# load manifest
mnf_df = loadManifestDf(manifest)
# get sample and index from cram filename
sample_id = cram_fl.split('/')[-1].split('#')[0]
idx = cram_fl.split('/')[-1].split('#')[-1].split(".cram")[0].replace('_','')
# find reference tag
ref_id = mnf_df.loc[(mnf_df["lims_id"]==sample_id) & (mnf_df['index']==int(idx))]["ref"].values[0]

# get fasta files for the reference tag
if ref_id == "PFA_GRC1_v1.0":
    fasta_wldcard = ref_database+"grc1/Pf_GRC1v1.0.fasta*"

if ref_id == "PFA_GRC2_v1.0":
    fasta_wldcard = ref_database+"grc2/Pf_GRC2v1.0.fasta*"

if ref_id == "PFA_Spec":
    fasta_wldcard = ref_database+"spec/Spec_v1.0.fasta*"

# copy to the working dir
system(f"ln -sf {fasta_wldcard} {getcwd()}")
# TODO work with symlinks instead of actually copy the files
