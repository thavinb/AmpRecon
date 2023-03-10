import json
import argparse

# --| FUNCTIONS |--------------------------------------------------------------
# load barcode definition

## Write COI INPUT functions

def loadBarcodeDef(barcode_json_path:str) -> dict:
    """
    Load barcode defintion from a json file

    ex:
    { "barcoding":{
        "barcode_ref": {
            "001": {"Chromosome": "Pf3D7_02_v3", "Locus": 376222, "RefAllele": "A"},
            "002": {"Chromosome": "Pf3D7_02_v3", "Locus": 470013, "RefAllele": "G"},
            "003": {"Chromosome": "Pf3D7_03_v3", "Locus": 656861, "RefAllele": "T"}
            }
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
    cur_keys = list(barcode_def_dct["barcoding"]["barcode_ref"].keys())
    print(barcode_def_dct["barcoding"]["barcode_ref"])
    for d_i in cur_keys:
        barcode_def_dct["barcoding"]["barcode_ref"][int(d_i)] = barcode_def_dct["barcoding"]["barcode_ref"].pop(d_i)

    return barcode_def_dct["barcoding"]['barcode_ref']

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
    four can be assigned:
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
        assert(char in ALLOWED_CHAR),f"ERROR: {char} is not a valid barcode character"
            
        
    coi_numbers = []
    for i, char_i in enumerate(barcode.upper()):
        pos_i = i+1
        # get barcode definition for pos_i
        ref_base_i = barcode_def_dct[pos_i]["RefAllele"]
        # sanity check
        assert(ref_base_i != "N"),"ERROR: reference allele set as N is forbidden. Sorry, it is the law!"
        
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

## Write coi GRC options

# Parser COI output to GRC
def parseCOIout(coiout_path:str)-> dict:
    """
    load McCOIL output summary file as dictionary

    Parameters
    ---

    coiout_path:<str>
        COI output summary file path
    
    coi output summary looks like:
    file\tCorP\tname\tmean\tmedian\tsd\tquantile0.025\tquantile0.975
    COIout\tC\tSPT59269\t1\t1\t0\t1\t1
    COIout\tC\tSPT42677\t2\t2\t0.48162\t1\t2
    COIout\tC\tSPT31022\t3\t3\t0.13761\t3\t3
    COIout\tC\tSPT32861\t2\t2\t0.00636\t2\t2
    """
    with open(coiout_path, "r") as coi_output:
        dct_out = {}
        is_header = True
        for line in coi_output:
            # load headers
            if is_header == True:
                col_names_lst = line.replace("\n","").split("\t")
                # be sure the columns that matters are on the expected places
                assert(col_names_lst[2] == "name")
                assert(col_names_lst[3] == "mean")
                is_header = False
                continue
            # get row values
            data_ln = line.replace("\n","").split("\t")
            # sanity check
            assert(len(data_ln) == len(col_names_lst))
            
            # skip sites P
            # we only want COI estimation for the "individuals"
            # so we only care about rows with CorP set as C
            if data_ln[1] == "P":
                continue
            
            # parse data to dict
            dct_out[data_ln[2]] = {col_names_lst[i]:data_ln[i] for i in range(1, len(col_names_lst))}
    return dct_out


# write GRC
def writeCOIgrc(coi_out_dct:dict, grc_out_flpath:str) -> None:
    """
    write GRC file from McCOIL output data

    Parameters
    ---
    coi_out_dct:<dict>
        A dictionary containing COI output summary data
        ex:{"SPT59269": {"CorP": "C", "name": "SPT59269", 
                         "mean": "1", "median": "1", 
                         "sd": "0", "quantile0.025": "1",
                          "quantile0.975": "1"},
            "SPT42677": {"CorP": "C", "name": "SPT42677", 
                        "mean": "2", "median": "2",
                        "sd": "0.48162","quantile0.025": "1",
                        "quantile0.975": "2"}
            }
    grc_out_flpath:<str>
        Path for the grc out file to be written
    """
    # get sample names
    samples_nms = list(coi_out_dct.keys())
    
    # write grc lines 
    with open(grc_out_flpath,"w") as coi_grc_fl:
        # header
        coi_grc_fl.write("ID\tMcCOIL\n")
        # rows
        for sample in samples_nms:
            coi_grc_fl.write(f"{sample}\t{coi_out_dct[sample]['mean']}\n")

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # get args
    parser = argparse.ArgumentParser(
        prog = "grc_process_mccoil_io.py",
        description = """
        A script to write McCOIL input from barcodes tsv files and write grc files from McCOIL output
        """
    )
    # Modes
    parser.add_argument("-write_mccoil_in", help="""
        prepare McCOIL input files from barcodes""", 
        action=argparse.BooleanOptionalAction)
    parser.add_argument("-write_coi_grc", help="""
        prepare McCOIL input files from barcodes""",
        action=argparse.BooleanOptionalAction)
    
    # inputs for write_mccoil_in mode
    parser.add_argument("--barcodes_files", nargs='+', help="""
        Path(s) to barcode tsv file for a batch of samples
    """, default=None)
    parser.add_argument("--config", help="""
        Path to a json file with barcodes definitions
    """, default=None)

    # inputs for write_coi_grc mode
    parser.add_argument("--mccoil_sum_file", nargs='+', help="""
        Path to McCOIL summary output file
    """, default=None)
    # output file name
    parser.add_argument("--output_file", default=None, help="""
        Path for the McCOIl input or grc file to be written 
        (default: <./McCOIL_in.tsv | ./coi.grc >).
        If more than one input file is provided, the output will have the same
        name of the respective input file added of a suffix ("_mccoil_in.tsv" | "_coi.grc")  
    """)
    
    args = {k:v for k,v in vars(parser.parse_args()).items()}
    
    print(" >> grc_process_mccoil_io.py <<")
    
    # handle modes
    try:
        assert((args["write_mccoil_in"]==None) or (args["write_coi_grc"]==None))
    except(AssertionError):
        print("ERROR: no mode of execution was specified. User must provide either '-write_mccoil_in' or '-write_coi_grc'.")
        print("       Choose wisely.")
        exit(1)
    try:
        assert((args["write_mccoil_in"]!= None) or (args["write_coi_grc"]!=None))
    except(AssertionError):
        print("ERROR: Only one mode of execution must be specified. You cannot have everything.")
        print("       User must provide either '-write_mccoil_in' or '-write_coi_grc'.")
        exit(1)


    # handle output default values
    if args["output_file"] == None:
        if args["write_mccoil_in"] == True:
            args["output_file"] = "./McCOIL_in.tsv"
        if args["write_coi_grc"] == True:
            args["output_file"] = "./coi.grc"

    if args["write_mccoil_in"]:
        # check if all required args were provided
        print(" -> Write McCOIL Input Mode On <-")
        if args["barcodes_files"] == None:
            print("    ERROR: 'barcodes_files' must be provided for this mode")
            exit(1)
        if args["barcode_def_file"] == None:
            print("     ERROR: 'barcodes_def_file' must be provided for this mode")
            exit(1)
        
        # if more than one, set output as same as input name
        total_batchs = len(args["barcodes_files"])
        print("    @ loading barcode definition...")
        barcode_def_dct = loadBarcodeDef(args["barcode_def_file"])
            
        for batch in args["barcodes_files"]:
            # load input data
            print(f"    @ loading samples barcodes ({batch})...")
            barcodes_samples_dct = loadSamplesBarcode(batch)
            print("    @ assigning SNP numbers for McCOIL...")
            assignSNPnumbers(barcodes_samples_dct)
            
            print("    @ writing McCOIL input tsv file...")
            if total_batchs == 1:
                writeMcCOILat(args["output_file"], barcode_def_dct, barcodes_samples_dct)
            
            if total_batchs > 1:
                # get batch simple file name (anything before the '.')
                prefix = batch.split("/")[-1].split(".")[0]
                # write COI input
                writeMcCOILat(prefix+'_mccoil_in.tsv', barcode_def_dct, barcodes_samples_dct)
        print(":: DONE ::")
        exit(0)

    if args["write_coi_grc"]:
        # check if all require args were provided
        print(" -> Write COI GRC Mode On <-")
        if args["mccoil_sum_file"] == None:
            print("    ERROR: No mccoil summary file was provided.")
            exit(1)

        total_batchs = len(args["mccoil_sum_file"])
        for batch in args["mccoil_sum_file"]:
            print(f"    @ loading {batch}")
            coi_out_dct = parseCOIout(batch)
            print(f"    @ writing {args['output_file']}...")
            if total_batchs == 1:
                writeCOIgrc(coi_out_dct, args["output_file"])
            if total_batchs > 1:
                # get batch simple file name (anything before the '.')
                prefix = batch.split("/")[-1].split(".")[0]
                # write COI grc
                writeCOIgrc(coi_out_dct, prefix+"_coi.grc")
        print(":: DONE ::")
        exit(0)