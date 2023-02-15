import csv
import os
import json

class AminoAcidCaller: #better name 

	def __init__(self, sample_id, genotypes_file, config_file):
		#need to think about how the process begins, does one file go in
		with open(config_file) as conf:
			self.config = json.load(conf)
		self.drl_info = self._read_drl_info()
		self.codon_key = self._read_codon_key()
		self.genotypes_files = self._read_genotypes_file(genotypes_file) #could probably be a dict with structure {sample_id: data???}
		self.sample_ids = self._get_sample_id(genotypes_file)
		self.complement_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
		self.haplotypes = {} #GRC1 grc2 

	@staticmethod
	def _read_tsv_row_generator(file):
		if not os.path.isfile(file):
			raise Exception(f"could not find {file}")
		with open(file) as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				yield row

	def _get_call(self, row):
		call = row['Gen']
		return call

	def _get_sample_id(self, file):
		sample_ids = set()
		with open(file) as f:
			for line in f:
				sample_ids.add(line[0])
		return sample_ids

	def _read_genotypes_file(self,file):
		'''
		a method that takes in a list of genotype files and parses them to a dictionary
		with {sample_id : {}}
		'''

		genotypes = { 'key' : {}, 'depth' : {} } #outer key 'key' is redundant

		for row in self._read_tsv_row_generator(file): #this is fine 
			chr = row['Chr']
			sample_id = row['Sample_ID']
			call = self._get_call(row)

			if sample_id not in genotypes['key']:
				genotypes['key'][sample_id] = {}
				genotypes['depth'][sample_id] = {}

			if chr not in genotypes['key'][sample_id]:
				genotypes['key'][sample_id][chr] = {}
				genotypes['depth'][sample_id][chr] = {}

			if row['Depth'] == '-':
				row['Depth'] = 0

			genotypes['key'][sample_id][chr][row['Loc']] = call
			genotypes['depth'][sample_id][chr][row['Loc']] = int(row['Depth'])
		
		return genotypes


	def _read_drl_info(self):
		"""Reads the DRL info file and returns a dict."""
		#this is fine - may be able to rearrange it in a way which is more palatable to reference - don't think that is super important
		ret = {'pos' : {}, 'core_order' : [], 'exp_order' : [], 'cons' : {}, 'aa_loc' : {}, 'strand' : {}, 'genes' : []}
		filename = "./amplicon_data/DRLinfo.txt" #change later 
		for row in self._read_tsv_row_generator(filename):
			if row['GeneName'] not in ret['pos']:
				ret['pos'][row['GeneName']] = {}
				ret['strand'][row['GeneName']] = {}
			if row['AA_Pos'] not in ret['pos'][row['GeneName']]:
				ret['pos'][row['GeneName']][row['AA_Pos']] = {'1' : row['pos1'], '2' : row['pos2'], '3' : row['pos3']}
				genename_aapos = f"{row['GeneName']}:{row['AA_Pos']}"
				ret['genes'].append(row['GeneName'])
				if row['Core'] == 'Y':
					ret['core_order'].append(genename_aapos)
				ret['exp_order'].append(genename_aapos)
				if row['Chr'] not in ret['cons']:
					ret['cons'][row['Chr']] = {}
				ret['cons'][row['Chr']][row['pos1']] = row['key1']
				ret['cons'][row['Chr']][row['pos2']] = row['key2']
				ret['cons'][row['Chr']][row['pos3']] = row['key3']
				ret['aa_loc'][genename_aapos] = row['Chr']
				ret['strand'][row['GeneName']][row['AA_Pos']] = row['Strand'] 
			ret['genes'] = list(set(ret['genes']))
		
		return ret

	def _read_codon_key(self):
		#think this has been superceeded by codon_table.py
		filename = "./amplicon_data/codonKey.txt"
		if not os.path.isfile(filename):
			raise Exception(f"could not find {file}")
		ret = {}
		with open(filename) as file:
			reader = csv.reader(file, delimiter='\t')
			for row in reader:
				if row[0] not in ret:
					ret[row[0]] = {}
				if row[1] not in ret[row[0]]:
					ret[row[0]][row[1]] = {}
				ret[row[0]][row[1]][row[2]] = row[4]
		
		return ret

	def is_missing(self, base):
		# perhaps a bit pointless
		return base == '-'

	def call_haplotypes(self):
		drl = self.drl_info
		core_genes = drl['core_order']
		sample_id = self.sample_id

		genotype_information = self.genotypes_files

		for gene_amino in drl['exp_order']:

			double_het = False

			gene, amino = gene_amino.split(":")

			if gene not in self.haplotypes:
				self.haplotypes[gene] = ''
			if gene_amino not in self.haplotypes:
				self.haplotypes[gene_amino] = ''

			if gene == 'PfCRT' and amino == '74':
				self.haplotypes[gene] += 'V'

			position_1 = drl['pos'][gene][amino]['1']
			position_2 = drl['pos'][gene][amino]['2']
			position_3 = drl['pos'][gene][amino]['3']
			chromosome = drl['aa_loc'][gene_amino]
			base_1 = '-'
			base_2 = '-'
			base_3 = '-'
			
			if position_1 not in genotype_information['key'][sample_id][chromosome]:
			#this is horrible, change it - this is to check if an amino position had been filtered out in a previous step, but 
			#is not clear what it is doing 
				if gene_amino in core_genes:
					self.haplotypes[gene] += '-'
				self.haplotypes[gene_amino] = '-'	
				continue

			if drl['cons'][chromosome][position_1] != '-':
				base_1 = drl['cons'][chromosome][position_1]
			if drl['cons'][chromosome][position_2] != '-':
				base_2 = drl['cons'][chromosome][position_2]
			if drl['cons'][chromosome][position_3] != '-':
				base_3 = drl['cons'][chromosome][position_3]

			if base_1 == '-':
				base_1 = genotype_information['key'][sample_id][chromosome][position_1] #how can i get a sample id here
			if base_2 == '-':
				base_2 = genotype_information['key'][sample_id][chromosome][position_2]
			if base_3 == '-':
				base_3 = genotype_information['key'][sample_id][chromosome][position_3]

			if drl['strand'][gene][amino] == '-':
				base_1 = self._complement(base_1)
				base_2 = self._complement(base_2)
				base_3 = self._complement(base_3)

			if self.is_missing(base_1) or self.is_missing(base_2) or self.is_missing(base_3):
				#this can go at some point
				self.haplotypes[gene_amino] = '-'
				if gene_amino in core_genes:
					self.haplotypes[gene] += '-'
				continue

			if len(base_1) == 1 and len(base_2) == 1 and len(base_3) == 1:
				#basic case 
				aa = self.translate_codon(base_1, base_2, base_3)
				if gene_amino in core_genes:
					self.haplotypes[gene] += aa
				self.haplotypes[gene_amino] = aa
				continue

			aa_1 = '-'
			aa_2 = '-'

			#the problem here is the double het cases cause the programme to crash, when checking if a single base is het yes it executes but we get a key error because we haven't
			#sufficiently dealt with the second het base

			#we should identify double het first

			if len(base_1) > 1 and len(base_2) > 1 or len(base_2) > 1 and len(base_3) > 1:
				double_het = True
				try:
					het_call = self._get_het_match_from_config(gene, amino, base_1, base_2, base_3)
					if gene_amino in core_genes:
						self.haplotypes[gene] += het_call
					self.haplotypes[gene_amino] += het_call
					continue
				except (TypeError, KeyError):
					if gene_amino in core_genes:
						self.haplotypes[gene] += '-'
					self.haplotypes[gene_amino] += '-'
					print(f"Double heterozygous case present in data not in list of special cases accounted for by the pipeline, please check input data at the locus {gene_amino}")

			elif len(base_1) > 1:
				allel_1 = base_1.split(',')[0]
				allel_2 = base_1.split(',')[1]
				aa_1 = self.translate_codon(allel_1, base_2, base_3)
				aa_2 = self.translate_codon(allel_2, base_2, base_3)
			elif len(base_2) > 1:
				allel_1 = base_2.split(',')[0]
				allel_2 = base_2.split(',')[1]
				aa_1 = self.translate_codon(base_1, allel_1, base_3)
				aa_2 = self.translate_codon(base_1, allel_2, base_3)
			elif len(base_3) > 1:
				allel_3 = base_3.split(',')[0]
				allel_4 = base_3.split(',')[1]
				aa_1 = self.translate_codon(base_1, base_2, allel_1)
				aa_2 = self.translate_codon(base_1, base_2, allel_2)

		if aa_1 != '-' and aa_2 != '-':
			if aa_1 == aa_2:
				aa_call = aa_1
			else:
				aa_call = f"[{aa_1}/{aa_2}]"

			self.haplotypes[gene_amino] = aa_call

			if gene_amino in core_genes:
				self.haplotypes[gene] += aa_call
		return self.haplotypes


	def _get_het_match_from_config(self, gene, amino, base_1, base_2, base_3):
		aa_1 = '-'
		aa_2 = '-'

		for case in self.config["DR_DOUBLE_HETEROZYGOUS_CASES"][gene]:
			
			if case[0] == amino and case[1][0] == base_1 and case[1][1] == base_2 and case[1][2] == base_3:
				aa_1 = case[2]
				aa_2 = case[3]

				if aa_1 != '-' and aa_2 != '-':
					if aa_1 == aa_2:
						aa_call = aa_1
						return aa_call
					else:
						return f"[{aa_1}/{aa_2}]"
				continue


	def translate_codon(self, base_1, base_2, base_3):
		amino_acid = '-'
		codon_translations = self.codon_key		
		amino_acid = codon_translations[base_1][base_2][base_3]
		return amino_acid
		

	def _complement(self, sequence):
		complement_sequence = []

		for nucleotide in sequence:
			if nucleotide in self.complement_bases:
				complement_sequence.append(self.complement_bases[nucleotide])
		return ''.join(complement_sequence)


class DoubleHeterozygousException(Exception):
	"""
	for when we encounter a double het that is not a special case 
	"""
	pass
