import vcf
import json
from collections import defaultdict

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
        match_threshold=0.95
    ):
        self.vcf_file = vcf.Reader(filename=vcf_file)
        self.sample = self.vcf_file.samples[0]
        self.min_maf = min_maf
        self.match_threshold = match_threshold
        self.chrom_key = self._read_json(chrom_key)
        self.species_ref = self._read_json(species_ref)
        self.alleles_depth_dict = self._convert_coords_store_alleles_depth()
        self._manipulate_maf()
        self.merged_alleles = self._merge_alleles()
        self.matched_loci = self._match_loci()
        self.species_label = self._get_species_tag()

    @staticmethod
    def _read_json(fp):
        return json.load(open(fp))

    def _convert_coords_store_alleles_depth(self):
        out = defaultdict(dict)
        for record in self.vcf_file:
            comboNum = self.chrom_key.get(record.CHROM).get(str(record.POS))
            alleles = [str(i) for i in [record.REF]+record.ALT if i]
            try:
                call = record.genotype(self.sample)
                called_alleles = list(set([alleles[i] for i in [int(j) for j in call["GT"].split("/")]]))
                DP = call["DP"]
            except IndexError:
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
        return 1-speciesDepth/totDepth

    def _manipulate_maf(self):
        """
        Iterate through each position in dict, 
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
        out = defaultdict(set)
        for pos, d_depth in self.alleles_depth_dict.items():
            for species in d_species_convert.values():
                out[pos].update(d_depth[species]["Allele"])
        return {k:sorted(list(v)) for k,v in out.items()}

    def _match_loci(self):
        first_item = list(self.species_ref.keys())[0]
        out = {s:0 for s in self.species_ref[first_item]}

        tot = 0
        for pos, d_alleles in self.species_ref.items():
            tot+=1
            sample_alleles = self.merged_alleles.get(pos,[])
            for species, allele in d_alleles.items():
                if allele in sample_alleles:
                    out[species]+=1
        out = {k:v/tot for k,v in out.items()}
        return out
    
    def _get_species_tag(self):
        label = ",".join([d_species_convert_match[k] for k,v in self.matched_loci.items() if v>=self.match_threshold])

        if not label:
            label = "-"

        return label
