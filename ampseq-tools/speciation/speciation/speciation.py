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
        chrom_key,
        species_ref,
        min_maf=0.01,
        match_threshold=0.95,
        min_total_depth=5,
        het_min_allele_depth=10
    ):
        self.vcf_file = vcf.Reader(filename=vcf_file)
        self.sample = self.vcf_file.samples[0]
        self.min_maf = min_maf
        self.match_threshold = match_threshold
        self.min_total_depth = min_total_depth
        self.het_min_allele_depth = het_min_allele_depth
        self.chrom_key = self._read_json(chrom_key)
        self.species_ref = self._read_json(species_ref)

        self.alleles_depth_dict = self._convert_coords_store_alleles_depth()
        self._manipulate_maf()
        self.merged_alleles = self._merge_alleles()
        self.matched_loci = self._match_loci()
        self.species_label = self._get_species_tag()

    @staticmethod
    def _read_json(fp):
        """
        read in json object and return dictionary
        """
        return json.load(open(fp))

    def _get_call(self, call, all_alleles):
        """
        Replicates het filtering logic applied in core pipeline.
        For a given call:
        1) get total depth
        2) if total depth > min total depth (5 from production)
            a) get unique genotype bases
            b) if call is het (alleles list length > 1)
                i) get allele depth for each allele
                ii) create alleles list from all allleles with AD greater than het_min_allele_depth (10 from production)
            c) if call not het do not modify alleles list
        3) if total depth less, empty alleles list and DP = 0
        4) return alleles and DP
        """
        DP = call.data.DP
        if DP >= self.min_total_depth:
            pattern = "/|\|"
            gt = set([int(i) for i in re.split(pattern, call["GT"]) if i])
            if isinstance(call["AD"],list):
                ad = [int(i) for i in call["AD"]]
            else:
                ad = [int(call["AD"])]
            d_ad = dict(zip(all_alleles, ad))
            alleles = [str(all_alleles[i]) for i in gt]
            if alleles:
                if len(alleles)>1:
                    alleles = [a for a in alleles if d_ad[a]>=self.het_min_allele_depth]
        else:
            alleles = []
            DP = 0
        return alleles, DP


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
        for record in self.vcf_file:
            comboNum = self.chrom_key.get(record.CHROM).get(str(record.POS))
            ref = record.REF
            alts = [i for i in record.ALT if i]
            all_alleles = [str(i) for i in [ref]+alts]
            try:
                assert not record.FILTER
                call = record.genotype(self.sample)
                called_alleles, DP = self._get_call(call, all_alleles)
            except (IndexError, AssertionError):
                called_alleles = []
                DP = 0
            
            species = d_species_convert.get(record.CHROM.split("_")[-1])
            out[comboNum][species] = {"Allele":called_alleles,"DP":DP}

        out_copy = out.copy()
        for pos, d_gt in out_copy.items():
            for species in d_species_convert.values():
                if species not in d_gt:
                    out[pos][species] = {"Allele":[],"DP":0}
        out_copy = None
        return out

    @staticmethod
    def _calculate_maf(speciesDepth, totDepth):
        """
        Calculate Pf/Pv MAF
        """
        return 1-speciesDepth/totDepth

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
        iterable = self.alleles_depth_dict.copy()
        
        for pos, d_gt in iterable.items():
            species=None
            totDepth=d_gt["Pf"]["DP"]+d_gt["Pv"]["DP"]
            if totDepth>0:
                if d_gt["Pf"]["DP"]>d_gt["Pv"]["DP"]:
                    species="Pf"
                    DP=d_gt["Pf"]["DP"]
                elif d_gt["Pv"]["DP"]>d_gt["Pf"]["DP"]:
                    species="Pv"
                    DP=d_gt["Pv"]["DP"]
            
            if species:
                maf = self._calculate_maf(DP, totDepth)
                if maf<self.min_maf:
                    if species=="Pf":
                        self.alleles_depth_dict[pos]["Pv"]["Allele"] = []
                    elif species=="Pv":
                        self.alleles_depth_dict[pos]["Pf"]["Allele"] = []
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
            for species in d_species_convert.values():
                out[pos].update(d_depth[species]["Allele"])
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
        for pos, d_alleles in self.species_ref.items():
            sample_alleles = self.merged_alleles.get(pos,[])
            if sample_alleles:
                tot+=1
            for species, allele in d_alleles.items():
                if allele in sample_alleles:
                    out[species]+=1
        out = {k:v/tot for k,v in out.items()}
        return out
    
    def _get_species_tag(self):
        """
        Join together species tags which have a frequency greater than match_threshold (0.95 in production)
        """
        label = "/".join([d_species_convert_match[k] for k,v in self.matched_loci.items() if v>=self.match_threshold])

        if not label:
            label = "-"

        return label
