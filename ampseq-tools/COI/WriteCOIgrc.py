import argparse

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
        coi_grc_fl.write("Sample\tMcCOIL\n")
        # rows
        for sample in samples_nms:
            coi_grc_fl.write(f"{sample}\t{coi_out_dct[sample]['mean']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = "WriteCOIGRC.py",
        description = "A script to write McCOIL COI values as a GRC"
    )

    parser.add_argument("-mccoil_sum_file", help="""
        Path to McCOIL summary output file
    """, required=True)
    parser.add_argument("-output_file", help="""
        Path for the grc containing the COI per sample to be written (default: ./coi.grc)
    """, default="./coi.grc")
    args = {k:v for k,v in vars(parser.parse_args()).items()}
    print(f"@ loading {args['mccoil_sum_file']}")
    coi_out_dct = parseCOIout(args['mccoil_sum_file'])
    print(f"@ writing {args['output_file']}...")
    writeCOIgrc(coi_out_dct,args["output_file"])
    print(":: DONE ::")