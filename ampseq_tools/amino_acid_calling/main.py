from AminoAcid import AminoAcidCaller
import argparse

def main():

	#rudimentary implementation to run methods only 
	#the pipeline should:
	#take a list of unlimited length containing genotype files and iteratively process each one OR
	#take directly from stdin - ie. cat g.txt
	# how do we get sample id? this will be a column from the genotypes.txt file, but we need to parse it somehow 
	# we also need to think about collating the output data into tsv at the end 
	

	parser = argparse.ArgumentParser()
	parser.add_argument('--genotypes_file', nargs='+', required=True)
	parser.add_argument('--species_config', required=True)
	args = parser.parse_args()

	if args.genotypes_file:
		for file in args.genotypes_file:
			output_data =[]
			caller = AminoAcidCaller(file, args.species_config)
			out = caller.call_haplotypes()
			output_data.append(out)
			write_out_grcs(output_data, caller.config['AMPLICON_PANEL_DESIGN'])


def write_out_grcs(haplotype_data:list, panel_design, extended=True):

	grc_1 = open("grc_1.tsv", 'a')
	grc_1.write("Sample_ID\t")

	#write header 
	for gene in panel_design["GENES"][0]:
		grc_1.write(gene+'\t')

	grc_1.write('\n')

	for entry in haplotype_data:
		sample_ids = list(entry.keys())
		for id in sample_ids:
			data = entry[id]
			grc_1.write(f"{id}\t{data['PfCRT']}\t{data['PfDHFR']}\t{data['PfDHPS']}\t{data['PfEXO']}\t{data['PfMDR1']}\t{data['PGB']}\n")


	if extended:
		grc_2 = open("grc_2.tsv", 'a')
		

if __name__ == "__main__":
	main()

	
