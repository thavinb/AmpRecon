import csv
import os

class AminoAcidCaller:

	def __init__(self, sample_id, genotypes_file):
		#with open(config_file) as conf:
		#	self.config_file = conf
		self.sample_id = sample_id
		self.drl_info = self._read_drl_info()
		self.codon_key = self._read_codon_key()
		self.genotypes_files = self._read_genotypes_file(genotypes_file) #could probably be a dict with structure {sample_id: data???}
		self.complement_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
		self.haplotypes = {}

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
		if row['Filt'] != 'PASS':  #Not expecting this to ever be executed??
			call = '-'
		return call

	def _read_genotypes_file(self,file):
		'''
		a method that takes in a list of genotype files and parses them to a dictionary
		with {sample_id : {}}
		'''
		genotypes = { 'key' : {}, 'depth' : {} }
		sample_id = self.sample_id

		for row in self._read_tsv_row_generator(file):
			chr = row['Chr']
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
		return base == '-'

	def call_haplotypes(self):
		drl = self.drl_info
		core_genes = drl['core_order']
		sample_id = self.sample_id

		genotype_information = self.genotypes_files

		for gene_amino in drl['exp_order']:

			double_het = false

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
				#this is horrible, change it
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
				base_3 == genotype_information['key'][sample_id][chromosome][position_3]

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
				aa = self.translate_codon(base_1, base_2, base_3)
				if gene_amino in core_genes:
					self.haplotypes[gene] += aa
				self.haplotypes[gene_amino] = aa
				continue

			aa_1 = '-'
			aa_2 = '-'
			
			if len(base_1) > 1:
				allel_1 = list(base_1)[0]
				allel_2 = list(base_1)[1]
				aa_1 = self.translate_codon(allel_1, base_2, base_3)
				aa_2 = self.translate_codon(allel_2, base_2, base_3)
			elif len(base_2) > 1:
				allel_1 = list(base_2)[0]
				allel_2 = list(base_2)[1]
				aa_1 = self.translate_codon(base_1, allel_1, base_3)
				aa_2 = self.translate_codon(base_1, allel_2, base_3)
			elif len(base_3) > 1:
				allel_1 = list(base_3)[0]
				allel_2 = list(base_3)[1]
				aa_1 = self.translate_codon(base_1, base_2, allel_1)
				aa_2 = self.translate_codon(base_1, base_2, allel_2)

			else:
				#do some double heterozygous logic here
				double_het = true
				bases = [base_1, base_2, base_3]

				for base in bases:
					if len(base) > 1:
						allel_1 = list(base)[0]
						allel_2 = list(base)[1]
						


			if aa_1 != '-' and aa_2 != '-':
				if aa_1 == aa_2:
					aa_call = aa_1
				else:
					aa_call = f"[{aa_1}/{aa_2}]"

			self.haplotypes[gene_amino] = aa_call

			if gene_amino in core_genes:
				#when i get a chance move this code to the start of the loop so i don't re write it all the time
				self.haplotypes[gene] += aa_call

			if 'PfCRT' in self.haplotypes and self.haplotypes['PfCRT'] == '-V---':
				self.haplotypes['PfCRT'] == '-----' 

	
		return self.haplotypes

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
