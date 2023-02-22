import argparse

# Parser COI output to GRC
def parseCOIout(coiout_path:str)-> dict:
    """
    load McCOIL output summary file as dictionary
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

    parser.add_argument("-mccoil_sum_fl", help="""
        Path to McCOIL summary output file
    """, required=True)
    parser.add_argument("-output_flnm", help="""
        Path for the grc containing the COI per sample to be written (default: ./coi.grc)
    """, default="./coi.grc")
    args = {k:v for k,v in vars(parser.parse_args()).items()}
    print(f"@ loading {args['mccoil_sum_fl']}")
    coi_out_dct = parseCOIout(args['mccoil_sum_fl'])
    print(f"@ writing {args['output_flnm']}...")
    writeCOIgrc(coi_out_dct,args["output_flnm"])
    print(":: DONE ::")