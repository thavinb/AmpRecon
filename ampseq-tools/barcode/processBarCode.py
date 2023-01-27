import pandas as pd


def __isGenotypeDFSane(gen_df):
    """
    do sanity checks on sample genotype dataframe
    """
    # are all required columns present?
    required_cols = ["Sanger_ID","Amplicon","Pos", 
                     "Chr","Loc","Gen","IsHet","Depth",
                     "Ale1_Dep","Ale2_Dep","Filt"]
    assert(required_cols in gen_df.columns)
    # PS: we should change the name of Sanger_id to something more generic
    # PS2: we may consider rename the BarcodeOrder columns to be in sync with
    #      the ones on the sample input (whatever it is).
    # check sanity of individual columns
    pass

class barcoder:
    '''
    A class to handle the generation of barcodes for sample amplicon data.

    ...

    Attributes
    ----------
    barcode_def_path : str
        barcode definition file path

    barcode_def_df : pd.DataFrame
        barcode definition dataframe

    barcode_length : int
        number of barcode positions defined
    
    sound : str
        the sound that the animal makes
    num_legs : int
        the number of legs the animal has (default 4)

    Methods
    -------
    says(sound=None)
        Prints the animals name and what sound it makes
    """
    '''

    def __init__(self, barcode_def_path, min_depth_treshold) -> None:
        # -- SANITY CHECK --
        # assert the input file exists
        # ------------------

        # set attributes of the barcode definition
        self.barcode_def_path = barcode_def_path
        self.barcode_def_df = pd.read_csv(barcode_def_path, sep="\t")
        assert(self.__checkSanityOfBarcodeDef() == True)
        self.barcode_length = len(self.barcode_def_df)
        self.min_depth_treshold = min_depth_treshold
        # load contigs headers name
        self.__barcode_sample_headers = self.__getSampleBarcodeHeaders()

    def __getSampleBarcodeHeaders(self) -> str:
        """Contig names are used as columns names at the sample barcode"""
        contig_names = self.barcdode_def_df["Chrom_ID"].values
        return "ssID\tMalGEN_ID\tStudy"+'\t'.join(contig_names)+'\n'

    def __checkSanityOfBarcodeDef(self) -> bool:
        """
        check is essential barcode assumptions are true
        """
        # -- are all required columns present?
        required_cols = ["Chrom_ID","VarStart","VarStop","Plex",
                        "VarType","Chromosome","Locus","RefAllele",
                        "VarPos","Concat","Class","Order","Present"]
        
        # PS: - My guess is that "ChromID", "Plex", "VarStart", "RefAllele",
        #     "Order" are the ones actually needed. If that is the case,
        #     lets keep the number of columns required at minimum.
        #
        #      - we use "Concat", but is just Chromosome+Locus (which is what 
        #     is present on the sample), should consider just to use the same name
        #     of the column used on the barcode
        #
        #      - We should rename Plex to something like "panel_name"
        # -- check "Order" columns
        # are all int?
        # are all unique values? (no barode position should be assigned twice)
        pass

    def loadBarcodeDataAsDict(self, sample_genotype_path, panel_name) -> dict:
        """
        write barcode file for a given sample genotype file

        Input
        -----
        sample_genotype : path
            sample genotype file

        out_barcode_file : path
            barcode output file name
        
        
        Return
        ------
        None
        """
        def __selectBarcodesAmplicons(gen_row):
            if gen_row["Amplicon"] in self.barcode_def_df["Chrom_ID"].values:
                return True
            else:
                return False
        
        # -- SANITY CHECK --
        # assert the sample_genotype file exists
        # ------------------

        # parse genotype file
        # if sample is genotype file
        genotype_df = pd.read_csv(sample_genotype_path, sep="\t")
        __isGenotypeDFSane(genotype_df)

        # get sample id
        # on the genotype file, one single ID is repeated for each position row
        # TODO: ADD ERROR MESSAGE
        assert(len(genotype_df["Sanger_ID"].unique()) == 1)
        sample_id = genotype_df["Sanger_ID"].iloc[0]

        # create dcitionary to hold the data
        dct_row = {"Sanger_ID":sample_id, # equivalent to ssID
                   # maybe we should not require or add the columns bellow
                   # at least not on this module
                   #"MalGen_ID":genotype_df["MalGen_ID"].values[0],
                   #"Study":sample_row["Study"].values[0], study is not present on the genotype, but somehow it is on the barcode
                }
        
        # -- LOAD SAMPLE DATA VALUES ---------------------------------
        # select only rows with contigs present in the barcode (optimization?)
        sample_valid_rows_df = genotype_df.loc[genotype_df.apply(__selectBarcodesAmplicons, axis=1)]
        
        # handle in case no row is valid (all bar positions should be set as "X")
        [myBeatifulErrorHandlingCode]

        # iterate over valid rows
        for i in range(0, len(sample_valid_rows_df)):
            sample_row = sample_valid_rows_df.iloc[i]
            sample_pos_i = sample_row["Pos"]
            sample_contig_i = sample_row["Chrom_ID"]
            sample_call_i = sample_row["Gen"]
            sample_depth_i = sample_row["Depth"]
            # -- Check if is a barcode position
            # set selection conditions
            valid_bar_pos = (self.barcode_def_df["VarPos"] == sample_pos_i)
            valid_contig_pos = (self.barcode_def_df["Chrom_ID"] == sample_contig_i)
            valid_pannel_pos = (self.barcode_def_df["Plex"] == panel_name)
            # get barcode row
            barcode_row = self.barcode_def_df.loc[(valid_contig_pos) and (valid_bar_pos) and valid_pannel_pos]
            
            # if barcode row is empty, current sample row is not a barcode position
            if len(barcode_row) == 0:
                continue
            else:
                assert(len(barcode_row)==1)
                # for some reason, the barcode columns are named as chromossome:locus,
                # which is stored at the barcode order as collate
                barcode_col_nm = barcode_row["Concat"]
                
                # check if sample barcode call depth is higher than treshold
                if sample_depth_i >= self.min_depth_treshold:
                    dct_row[barcode_col_nm] = sample_call_i
                # if not, set barcode position as "X"
                else:
                    dct_row[barcode_col_nm] = "X"
        
        # add not valid positions
        # get Amplicon (sample) not available at ChromID (barcode)
        amplicons_at_sample = genotype_df["Amplicon"].unique()
        
        not_in_sample_chromID = [x for x in amplicons_at_sample if x not in self.barcode_def_df["chromID"].values]
        for empty_chromID in not_in_sample_chromID:
            barcode_rows = self.barcode_def_df.loc[self.barcode_def_df["chromID"] == empty_chromID]
            for i in range(0,len(barcode_rows)):
                pos_row = barcode_rows.iloc[i]
                barcode_col_nm = pos_row["Concat"]
                dct_row[barcode_col_nm] = "X"

        return dct_row

    def writeBarcodesForSamplesBatch(self, sample_lst, output_flpth) -> pd.DataFrame:
        """
        Given a list os sample files (genotyping), write
        """
        ## WARN WE NEED TO ADD THE PANNEL NAME TO GENOTYPE FILE
        
        rows_dct_lst = []
        for sample_flpath in sample_lst:
            rows_dct_lst.append(self.loadBarcodeDataAsDict(sample_flpath))
        
        barcodes_df = pd.DataFrame(rows_dct_lst)
        barcodes_df.to_csv(output_flpth, sep="\t")
        return barcodes_df

if __name__ == '__main__':
    
    # arg parser
    # args: barcode_order_fl, min_depth_treshold, output_flpth, samples_in_fl_lst
    #       config_fl
    
    # if config_fl not None:
        # parse arguments from config_fl
    # else
        # get arguments from command line
    
    #barcode_def_in = "./barcodeOrder_vivax.tsv"
    #min_depth_treshold = 100 # need to check which one is used at curr production pipline and set as the default

    # create barcoder
    barcoder_obj = barcoder(barcode_def_path=barcode_def_in)
    # write barcodes file
    barcoder_obj.writeBarcodesForSamplesBatch(samples_in_fl_lst,output_flpth)
    


# get barcode definition tsv
#barcode_def_in = "./barcodeOrder.tsv"
#barcode_def_in = "./barcodeOrder_vivax.tsv"

#sample_genotype_fl = "./4977STDY8014273.tsv"
