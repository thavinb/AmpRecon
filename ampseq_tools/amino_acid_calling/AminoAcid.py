
class AminoAcidCaller:

	def __init__(self, sample_id, config_file, genotypes_file):
		with open(config_file) as conf:
			self.config_file = conf

		self.drl_info = _read_drl_info()
		self.codon_key = _read_codon_key()
		self.genotypes_file = _read_genotypes_file()
		self.haplotypes = {}

	def _get_core_haplotype():
		'''
		a method given some input paramters, generates a list of genes that make 
		up the core haplotype positions
		'''
		pass

	def _get_ref_alleles(gene, amino_position):
		'''
		a method given a gene and amino position, generates the nucleotide position
		numbers for each base in the codon
		'''
		pass

	def _do_basecall_replacement():
		'''
		a method given some input parameters, will replace empty calls in the DRL with
		experimentally derived calls
		'''
		pass

	def _is_reverse_stranded(gene):
		'''
		a method given a gene, will return true if the gene is on a reverse strand
		'''
		pass

	def _impute_pfCRT73():
		'''
		a method that will append a V to the haplotypes dict 
		'''
		pass

	def _is_missing(base):
		'''
		a method that will return true if a base equals '-'
		'''
		pass

	def _translate_codon(base1, base2, base3): 
		'''
		a method that will look up the codon information and return its constituent
		amino acid single letter code
		'''
		pass

	def _get_aa_match():
		'''
		a method that will look up whether the data concerns itself with a particular 
		het case in matches_bas conf object
		'''
		pass

	def write_haplotype():
		'''
		a method that will write a dictionary to tsv, the format should be
		sample_id\tgene\thaplotype
		'''
		pass

