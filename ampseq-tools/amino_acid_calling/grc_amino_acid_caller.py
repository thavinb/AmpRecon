import csv
import os
import json
import argparse


class AminoAcidCaller: #better name 

	def __init__(self, genotypes_files, config_file, drl_information_file, codon_key):
		#need to think about how the process begins, does one file go in
		with open(config_file) as conf:
			self.config = json.load(conf)
		self.drl_info = self._read_drl_info(drl_information_file)
		self.codon_key = self._read_codon_key(codon_key)
		self.genotypes_files = genotypes_files
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

	def _read_genotypes_file(self, files: list):
		'''
		a method that takes in a list of genotype files and parses them to a dictionary
		with {sample_id : {}}
		'''
		genotypes = {}

		for genotype_file in files:
			for row in self._read_tsv_row_generator(genotype_file):  
				chr = row['Chr']
				sample_id = row['MalGen_ID']
				call = row['Gen'] 

				if sample_id not in genotypes:
					genotypes[sample_id] = {}

				if chr not in genotypes[sample_id]:
					genotypes[sample_id][chr] = {}

				genotypes[sample_id][chr][row['Loc']] = call
		
		return genotypes


	def _read_drl_info(self, filename):
		"""Reads the DRL info file and returns a dict."""
		#this is fine - may be able to rearrange it in a way which is more palatable to reference - don't think that is super important
		ret = {'pos' : {}, 'core_order' : [], 'exp_order' : [], 'cons' : {}, 'aa_loc' : {}, 'strand' : {}, 'genes' : []}
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
			ret['genes'] = list(dict.fromkeys(ret['genes']))
		return ret

	def _read_codon_key(self, filename):
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

		genotype_information = self._read_genotypes_file(self.genotypes_files) 
		sample_id = list(genotype_information.keys())

		for id in sample_id:
			self.haplotypes[id] = {}
			grc_1_cols = {}
			grc_2_cols = {}

			for gene_amino in drl['exp_order']:
				gene, amino = gene_amino.split(":")
				#begin by creating empty dictionary entries for the gene (grc1)	and gene_amino (grc2)			
				if gene not in self.haplotypes[id]:
					self.haplotypes[id][gene] = []
					grc_1_cols[gene] = True
				
				if gene in self.config['grc_amino_acid_caller']['DR_HAPLOTYPE_SPECIAL_CASES']:
					if amino in self.config['grc_amino_acid_caller']['DR_HAPLOTYPE_SPECIAL_CASES'][gene]:
						aa = self.config['grc_amino_acid_caller']['DR_HAPLOTYPE_SPECIAL_CASES'][gene][amino]
						if gene_amino in core_genes:
							self.haplotypes[id][gene].append(aa)
						continue


				if gene_amino not in self.haplotypes[id]:
					self.haplotypes[id][gene_amino] = []
					grc_2_cols[gene_amino] = True

				#get nucleotide positions from drlinfo.txt
				position_1 = drl['pos'][gene][amino]['1']
				position_2 = drl['pos'][gene][amino]['2']
				position_3 = drl['pos'][gene][amino]['3'] 
				chromosome = drl['aa_loc'][gene_amino]

				base_1 = '-'
				base_2 = '-'
				base_3 = '-'

				#check for missing position in genotypes file we do position one as a shorthand - if its missing it doesn't matter
				#if the other bases are present or not
				if chromosome not in genotype_information[id]:
					if gene_amino in core_genes:
						self.haplotypes[id][gene].append('-')
					self.haplotypes[id][gene_amino].append('-')
					continue

				#assign reference bases to positions, if they exist in the drlinfo.txt
				if drl['cons'][chromosome][position_1] != '-':
					base_1 = drl['cons'][chromosome][position_1]
				if drl['cons'][chromosome][position_2] != '-':
					base_2 = drl['cons'][chromosome][position_2]
				if drl['cons'][chromosome][position_3] != '-':
					base_3 = drl['cons'][chromosome][position_3]

				if base_1 == '-' and position_1 in genotype_information[id][chromosome]:
					base_1 = genotype_information[id][chromosome][position_1]
					if drl['strand'][gene][amino] == '-':
						base_1 = self._complement(base_1)
				if base_2 == '-' and position_2 in genotype_information[id][chromosome]:
					base_2 = genotype_information[id][chromosome][position_2]
					if drl['strand'][gene][amino] == '-':
						base_2 = self._complement(base_2)
				if base_3 == '-' and position_3 in genotype_information[id][chromosome]:
					base_3 = genotype_information[id][chromosome][position_3]
					if drl['strand'][gene][amino] == '-':
						base_3 = self._complement(base_3)

				#if the data is on the reverse strand, we change the base to its complement

				#if any of the positions are missing, we call a missing aa and move on
				if self.is_missing(base_1) or self.is_missing(base_2) or self.is_missing(base_3):
					if gene_amino in core_genes:
						self.haplotypes[id][gene].append('-')
					self.haplotypes[id][gene_amino].append('-')
					continue

				#basic case of translating three homozygous bases
				if len(base_1) == 1 and len(base_2) == 1 and len(base_3) == 1: 
					aa = self.translate_codon(base_1, base_2, base_3)
					if gene_amino in core_genes:
						self.haplotypes[id][gene].append(aa)
					self.haplotypes[id][gene_amino].append(aa)
					continue

				#if we come here we must have a heterozguous

				aa_1 = '-'
				aa_2 = '-'

				#do the double het first - and then do the single het 

				#we should identify double het first

				if len(base_1) > 1 and len(base_2) > 1 or len(base_1) > 1 and len(base_3) > 1 or len(base_2) > 1 and len(base_3) > 1:
					try:
						het_call = self._get_het_match_from_config(gene, amino, base_1, base_2, base_3)
						assert(het_call != None)
						if gene_amino in core_genes:
							self.haplotypes[id][gene].append(het_call)
						self.haplotypes[id][gene_amino].append(het_call)
						continue

					except(AssertionError):
						if gene_amino in core_genes:
							self.haplotypes[id][gene].append('-')
						self.haplotypes[id][gene_amino].append('-')
						print(f"Double heterozygous case present in data not in list of special cases accounted for by the pipeline, please check input data at the locus: {gene_amino} for sample: {id}")
						continue

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
					allel_1 = base_3.split(',')[0]
					allel_2 = base_3.split(',')[1]
					aa_1 = self.translate_codon(base_1, base_2, allel_1)
					aa_2 = self.translate_codon(base_1, base_2, allel_2)
				else:
					raise HaplotypeProcessingException

				if aa_1 != '-' and aa_2 != '-':
					if aa_1 == aa_2:
						aa_call = aa_1
					else:
						aa_call = f"[{aa_1}/{aa_2}]"
					self.haplotypes[id][gene_amino].append(aa_call)
					if gene_amino in core_genes:
						self.haplotypes[id][gene].append(aa_call)

		for sample in self.haplotypes:
			for gene_amino in drl['exp_order']:
				gene, amino = gene_amino.split(":")
				hap = self.haplotypes[sample][gene]
				if gene in self.config['grc_amino_acid_caller']['DR_HAPLOTYPE_SPECIAL_CASES']:
					num_hc_for_gene_in_conf = len(self.config['grc_amino_acid_caller']['DR_HAPLOTYPE_SPECIAL_CASES'][gene])
					num_non_blank_for_gene = 0
					for call in self.haplotypes[sample][gene]:
						if call != '-':
							num_non_blank_for_gene += 1
					if num_hc_for_gene_in_conf == num_non_blank_for_gene:
						for i in range(len(self.haplotypes[sample][gene])):
							self.haplotypes[sample][gene][i] = '-'


		output_data = {

				"grc_1_cols" : grc_1_cols,
				"grc_2_cols" : grc_2_cols,
				"data" : self.haplotypes

		}
		return output_data

	def _get_het_match_from_config(self, gene, amino, base_1, base_2, base_3):
		aa_1 = '-'
		aa_2 = '-'

		if gene in self.config['grc_amino_acid_caller']['DR_DOUBLE_HETEROZYGOUS_CASES']:
			for case in self.config['grc_amino_acid_caller']["DR_DOUBLE_HETEROZYGOUS_CASES"][gene]:
				
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
			if nucleotide == '-':
					complement_sequence.append('-')
					break
			if nucleotide in self.complement_bases:
				complement_sequence.append(self.complement_bases[nucleotide])

		if len(complement_sequence) > 1:
			return ','.join(complement_sequence)
			
		return ''.join(complement_sequence)


class HaplotypeProcessingException(Exception):
	def __init__(self, message= "Error in processing genotype information, please check the input data"):
		self.message = message
		super.__init__(self.message)

class InputAmpliconFormattingException(Exception):
	def _init__(self, message= "Data contains nucleotide information that was not expected, expected base data is: A,T,G,C or '-'"):
		self.message = message
		super._init__(self.message)


def write_out_grcs(haplotype_data:list, drug_resistance_loci, output_grc1_file, output_grc2_file, extended=True):

	grc_1 = open(f"{output_grc1_file}", 'w+')
	grc_1.write("ID\t")
	sample_id = list(haplotype_data['data'].keys())

	if extended:
		grc_2 = open(f"{output_grc2_file}", 'w+')
		grc_2.write("ID\t")
		positions = haplotype_data['grc_2_cols']
		
		grc_2.write('\t'.join(positions) + '\n')

		for id in sample_id:
			line = []
			data = haplotype_data["data"][id]
			line.append(id)
			for pos in positions:
				line.append(data[pos][0])
			grc_2.write('\t'.join(line) + '\n')
		

	grc_1.write('\t'.join(drug_resistance_loci['genes']))
	grc_1.write('\n')

	for id in sample_id:
		line = []
		data = haplotype_data["data"][id]
		line.append(id)
		for gene in drug_resistance_loci['genes']:
			line.append(data[gene])

		dat = []
		for hap in line:
			hap = ''.join(hap)
			dat.append(hap)
		grc_1.write('\t'.join(dat) + '\n')
	
if __name__ == "__main__":
	#rudimentary implementation to run methods only 
	#the pipeline should:
	#take a list of unlimited length containing genotype files and iteratively process each one OR
	#take directly from stdin - ie. cat g.txt
	# how do we get sample id? this will be a column from the genotypes.txt file, but we need to parse it somehow 
	# we also need to think about collating the output data into tsv at the end 
	

	parser = argparse.ArgumentParser()
	parser.add_argument('--genotype_files', nargs='+', required=True)
	parser.add_argument('--config', required=True)
	parser.add_argument('--output_grc1_file', required=True) #add grc2 parameter 
	parser.add_argument('--output_grc2_file', required=True)
	parser.add_argument('--drl_information_file', required=True)
	parser.add_argument('--codon_key_file', required=True)

	args = parser.parse_args()

	grc_1 = open(f"{args.output_grc1_file}", 'w+')
	grc_2 = open(f"{args.output_grc2_file}", 'w+')
	caller = AminoAcidCaller(args.genotype_files, args.config, args.drl_information_file, args.codon_key_file)
	output = caller.call_haplotypes()
	write_out_grcs(output, caller.drl_info, args.output_grc1_file, args.output_grc2_file)

