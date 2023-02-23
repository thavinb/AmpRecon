from collections import defaultdict

d_species_convert = {
    "falciparum":"Pf",
    "vivax":"Pv"
}

class Speciate:
    def __init__(
        self,
        d_genotype_file: dict,
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
