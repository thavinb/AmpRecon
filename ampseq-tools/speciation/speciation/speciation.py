import vcf
import json
from collections import defaultdict
import re

d_species_convert = {
    "falciparum":"Pf",
    "vivax":"Pv"
}

d_species_convert_match = {
    "Falciparum":"Pf",
    "Vivax":"Pv",
    "Knowlesi":"Pk",
    "Ovale":"Po",
    "Malariae":"Pm"
}

class Speciate:
    def __init__(
        self,
        vcf_file,
        barcode,
        chrom_key,
        species_ref,
        min_maf=0.01,
        match_threshold=0.95,
        min_total_depth=10,
        het_min_allele_depth=5,
        default_species="Pf",
        default_species_barcode_coverage=51,
        d_species_convert=d_species_convert,
        d_species_convert_match=d_species_convert_match
    ):
        #initialise variables
        self.vcf_file = vcf.Reader(filename=vcf_file)
        self.sample = self.vcf_file.samples[0]
        self.min_maf = min_maf
        self.match_threshold = match_threshold
        self.min_total_depth = min_total_depth
        self.het_min_allele_depth = het_min_allele_depth
        self.chrom_key = self._read_json(chrom_key)
        self.species_ref = self._read_json(species_ref)
        self.d_species_convert = d_species_convert
        self.d_species_convert_match = d_species_convert_match
        self.write_out = {}

        #run main logic flow
        self.alleles_depth_dict = self._convert_coords_store_alleles_depth()
        self._manipulate_maf()
        self.merged_alleles = self._merge_alleles()
        self.matched_loci = self._match_loci()
        self.species_label = self._get_species_tag(
            barcode, 
            default_species, 
            default_species_barcode_coverage
            )

    @staticmethod
    def _read_json(fp):
        """
        read in json object and return dictionary
        """
        if isinstance(fp, str):
            return json.load(open(fp))
        elif isinstance(fp, dict):
            return fp

    def _get_call(self, call, all_alleles):
        """
        Replicates het filtering logic applied in core pipeline.
        For a given call:
        1) get total depth
        2) if total depth > min total depth (10 from production)
            a) get unique genotype bases
            b) if call is het (alleles list length > 1)
                i) get allele depth for each allele
                ii) create alleles list from all allleles with AD greater than het_min_allele_depth (5 from production)
            c) if call not het do not modify alleles list
        3) if total depth less, empty alleles list and DP = 0
        4) return alleles and DP
        """
        DP = call.data.DP
        #check total depth above minimum threshold
        if DP >= self.min_total_depth:
            #split genotype irrespective of phasing, get genotyped positions
            pattern = "/|\|"
            gt = set([int(i) for i in re.split(pattern, call["GT"]) if i])

            #place ad into a list whether is is a single value or not
            if isinstance(call["AD"],list):
                ad = [int(i) for i in call["AD"]]
            else:
                ad = [int(call["AD"])]
            #create allele -> ad dictionary
            d_ad = dict(zip(all_alleles, ad))
            #create list of all genotyped alleles - indexes from gt
            alleles = [str(all_alleles[i]) for i in gt]

            #If allele not genotyped, remove from ad list, remove ad from total depth
            for a in all_alleles:
                if a not in alleles:
                    ad1 = d_ad[a]
                    ad.remove(ad1)
                    DP-=ad1

            alleles_out = alleles.copy()
            for a in alleles:
                ad1 = d_ad[a]
                #if allele has depth below het_min_allele_depth, remove from alleles_out, remove AD for allele,
                #subtract AD from DP
                if ad1<self.het_min_allele_depth:
                    alleles_out.remove(a)
                    ad.remove(ad1)
                    DP-=ad1
        #Need to apply this check explicitely again here in case DP is 
        #now below threshold after above manipulation
        if not DP >= self.min_total_depth:
            alleles_out = []
            DP = 0
        return alleles_out, DP


    def _convert_coords_store_alleles_depth(self):
        """
        For each record in spec vcf file
        1) get comboNum coord from chrom_key dict
        2) ensure that record not filtered and that it has a genotype for the sample
           (otherwise give empty list and DP = 0)
        3. get called alleles and DP from _get_call
        4. get the species name from the contig identifier, transalate with d_species_convert
        5. add the called alleles and depth to the out dict
        6. Set empty entries where either genotype is missing
        """
        out = defaultdict(dict)
        #iterate through VCF object
        for record in self.vcf_file:
            #get "comboNUM" position from chrom_key
            comboNum = self.chrom_key.get(record.CHROM,{}).get(str(record.POS))

            #get ref and alts and combine in to all_alleles list
            ref = record.REF
            alts = [i for i in record.ALT if i]
            all_alleles = [str(i) for i in [ref]+alts]
            try:
                #get the genotype call for the record and send to get call method
                call = record.genotype(self.sample)
                called_alleles, DP = self._get_call(call, all_alleles)
            except IndexError as e:
                #if IndexError raised then sample missing genotype call at that position
                called_alleles = []
                DP = 0
            
            #get the aligned species simple name (Pf/Pv) and write to out dictionary
            species = self.d_species_convert.get(record.CHROM.split("_")[-1])
            out[comboNum][species] = {"Allele":called_alleles,"DP":DP}

        #Iterate through and apply empty rows to any position missing a species call
        out_copy = out.copy()
        for pos, d_gt in out_copy.items():
            for species in self.d_species_convert.values():
                if species not in d_gt:
                    out[pos][species] = {"Allele":[],"DP":0}
        out_copy = None
        return out

    @staticmethod
    def _calculate_maf(speciesDepth, totDepth):
        """
        Calculate Pf/Pv MAF
        """
        return 1-int(speciesDepth)/int(totDepth)

    def _manipulate_maf(self):
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
            if totDepth>0:
                if d_gt["Pf"]["DP"]>d_gt["Pv"]["DP"]:
                    species="Pf"
                    DP=d_gt["Pf"]["DP"]
                elif d_gt["Pv"]["DP"]>d_gt["Pf"]["DP"]:
                    species="Pv"
                    DP=d_gt["Pv"]["DP"]
            
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

    def _merge_alleles(self):
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

    def _match_loci(self):
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
    
    def _get_species_tag(self, barcode, default_species, min_coverage):
        """
        Join together species tags which have a frequency greater than match_threshold (0.95 in production).

        If barcode has sufficient coverage (51 non missing positions in production) append 
        default species (Pf in production)
        """
        #obtain list of all species with matched loci greater than threshold
        all_species_out = [self.d_species_convert_match[k] for k,v in self.matched_loci.items() if v>=self.match_threshold]

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
