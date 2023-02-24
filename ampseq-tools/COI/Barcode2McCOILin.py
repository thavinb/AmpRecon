import json
import argparse

# --| FUNCTIONS |--------------------------------------------------------------
# load barcode definition 
def loadBarcodeDef(barcode_json_path:str) -> dict:
    """
    Load barcode defintion from a json file

    ex:
    {
    "barcode_ref": {
        "001": {"Chromosome": "Pf3D7_02_v3", "Locus": 376222, "RefAllele": "A"},
        "002": {"Chromosome": "Pf3D7_02_v3", "Locus": 470013, "RefAllele": "G"},
        "003": {"Chromosome": "Pf3D7_03_v3", "Locus": 656861, "RefAllele": "T"}
        }
    }

    Parameters
    ---
    barcode_json_path:<path>
        path to a json file containing barcode definitions
    
    Returns
    ---
    dictionary with the barcode definitions
    """
    barcode_def_dct = json.load(open(barcode_json_path,"r"))

    # intify keys
    cur_keys = list(barcode_def_dct["barcode_ref"].keys())
    for d_i in cur_keys:
        barcode_def_dct["barcode_ref"][int(d_i)] = barcode_def_dct["barcode_ref"].pop(d_i)

    return barcode_def_dct['barcode_ref']


# load samples barcode file
def loadSamplesBarcode(samples_brcd_path:str)-> dict:
    """
    load tsv file with barcodes per sample

    ex:
    '''
    Sample\tBarcode
    PT59269\tAGTCAGTGCAAACGTCTAA
    SPT4416\tXXXXAXNXXAXXXXXXXXX
    '''
    """
    barcodes_samples_dct = {}
    with open(samples_brcd_path,"r") as barcode_infl:
        skip = True
        for line in barcode_infl:
            # skip header
            if skip == True:
                skip = False
                continue
            data_ln = line.replace("\n","").split("\t")
            # sanity check
            err_m = f"ERROR: more than 2 cols for sample {data_ln[0]} at {samples_brcd_path}"
            assert(len(data_ln)==2),err_m
        
            # {sample_id:{barcode:"NNNN"}}
            barcodes_samples_dct[data_ln[0]] = {"barcode":data_ln[1]}
    
    return barcodes_samples_dct

# assign SNP numbers for McCoil
def __getSNPnumbersFrom(barcode:str, barcode_def_dct:dict)-> list:
    """
    For The Real McCOIL categorical method, for each site analysed
    three four can be assigned:
        * -1 :if base is missing on a given individual/sample
        * 0.5: if base is heterozygous
        * 1 : if homozygous major allele
        * 0 : if homozygous minor allele
     
    This functions assumes the base encoding bellow:
        * "X", if position is missing
        * "N", if position is heterozygous (0.5)
        * same as reference, if position is homozygous minor
        * not same as reference, if position is homozygous major
    
    only A,T,C,G,X and N are accepted as valid bases
    
    Parameters
    ---

    barcode:<str>
        a given barcode string
    
    barcode_def_dct:dict
        barcode definition dictionary [{pos_i:{"RefAllele":refbase}}]
    
    Return
    ---
    <list>
        a list of the assigned numbers for each barcode position
        (same order of appearance)
    """
    # sanity check barcode
    ALLOWED_CHAR = ["X","N","A","T","C","G"]
    unique_char = list(set(barcode.upper())) 
    for char in unique_char:
        if char not in ALLOWED_CHAR:
            raise 
    coi_numbers = []
    for i, char_i in enumerate(barcode.upper()):
        pos_i = i+1
        # get barcode definition for pos_i
        ref_base_i = barcode_def_dct[pos_i]["RefAllele"]
        # sanity check
        assert(ref_base_i != "N"),"ERROR: reference allele set as N is forbiden. Sorry, is the law"
        
        # if "X" then is considered as a missing postitions, set -1
        # default assumption is assuming base is missing
        coi_i = -1
        # if "N" then is considered heterozygoes, set 0.5
        if char_i == "N":
            coi_i = 0.5
        # if equal to ref is considered homozygous minor,set 0
        if (char_i == ref_base_i) and (char_i != "X") and (char_i != "N"):
            coi_i = 0
        # if not equal to ref is considered homozygous major, set 1
        if (char_i != ref_base_i) and (char_i != "X") and (char_i != "N"):
            coi_i = 1
        coi_numbers.append(coi_i)

    return coi_numbers

def assignSNPnumbers(barcodes_samples_dct:dict) -> None:
    for sample_i in barcodes_samples_dct.keys():
        barcode_i = barcodes_samples_dct[sample_i]["barcode"]
        coi_in_i = __getSNPnumbersFrom(barcode_i, barcode_def_dct)
        barcodes_samples_dct[sample_i]["coi_in"] = coi_in_i

# write RealMcCoil input
def writeMcCOILat(out_flpath:str, barcode_def_dct:dict, barcodes_samples_dct:dict) -> None:
    """
    wirte COI input file from barcode definition dictionary

    Parameter:
    ---

    out_flpath:<str>
        output filepath

    barcode_def_dct:<dict>
        barcode definition dictionary
        ex:
        {
            1: {"Chromosome": "Pf3D7_02_v3", "Locus": 376222},
            2: {"Chromosome": "Pf3D7_02_v3", "Locus": 470013},
            2: {"Chromosome": "Pf3D7_03_v3", "Locus": 656861}
            }
        }

    barcodes_samples_dct:<dict>
        dictionary containing the barcodes 
        ex:
        {
        'PT59269':{'coi_in:'[0,0,0]},
        'SPT4416':{'coi_in':[-1,-1,-1]}
        'SPT4417':{'coi_in':[0,0,1]}
        }
    """
    with open(out_flpath, "w") as incoi_fl:
        # write col names
        header = "sample_id\t"
        for pos_i in barcode_def_dct.keys():
            contig_nm_i = barcode_def_dct[pos_i]["Chromosome"]
            locus_nm_i = str(barcode_def_dct[pos_i]["Locus"])
            site_header_i = contig_nm_i+"_"+locus_nm_i
            header+=site_header_i+"\t"
        incoi_fl.write(header[0:-1]+"\n")

        # write rows
        for sample in barcodes_samples_dct.keys():
            values = '\t'.join([str(x) for x in barcodes_samples_dct[sample]['coi_in']])
            incoi_fl.write(str(sample)+"\t"+values+"\n")

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # get args
    parser = argparse.ArgumentParser(
        prog = "Barcode2McCOILin.py",
        description = "A script to write McCOIL input from barcodes tsv files"
    )

    parser.add_argument("-barcodes_in", help="""
        Path to input barcode tsv file
    """, required=True)
    parser.add_argument("-barcode_def", help="""
        Path to a json file with barcodes definitions
    """, required=True)
    parser.add_argument("-output_flnm", default="./McCOIL_in.tsv", help="""
        Path for the McCOIl input file to be written (default: ./McCOIL_in.tsv)
    """)
    
    args = {k:v for k,v in vars(parser.parse_args()).items()}
    print(" >> Barcode2McCOILin.py <<")
    # load input data
    print("@ loading barcode definition...")
    barcode_def_dct = loadBarcodeDef(args["barcode_def"])
    print("@ loading samples barcodes...")
    barcodes_samples_dct = loadSamplesBarcode(args["barcodes_in"])
    print("@ assigning SNP numbers for McCOIL...")
    assignSNPnumbers(barcodes_samples_dct)

    print("@ writing McCOIL input tsv file...")
    writeMcCOILat(args["output_flnm"], barcode_def_dct,barcodes_samples_dct)
    print(":: DONE ::")
