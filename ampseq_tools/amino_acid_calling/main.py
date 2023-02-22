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
	parser.add_argument('--output_file', required=True)

	args = parser.parse_args()

	if args.genotypes_file:
		output_data =[]
		for file in args.genotypes_file:
			caller = AminoAcidCaller(file, args.species_config)
			out = caller.call_haplotypes()
			output_data.append(out)
			write_out_grcs(output_data, caller.drl_info, args.output_file)


def write_out_grcs(haplotype_data:list, drug_resistance_loci, output_file, extended=True):

	grc_1 = open(f"{output_file}", 'w')
	grc_1.write("ID\t")

	if extended:
		grc_2 = open(f"{output_file}.extended", 'w')
		grc_2.write("ID\t")
		positions = drug_resistance_loci['exp_order']
		grc_2.write('\t'.join(positions) + '\n')
		
		for entry in haplotype_data:
			line = []
			sample_ids = list(entry.keys())
			for id in sample_ids:
				data = entry[id]
				line.append(id)
				for pos in positions:
					line.append(data[pos])
				grc_2.write('\t'.join(line) + '\n')

	for gene in drug_resistance_loci['genes']:
		grc_1.write(gene+'\t')

	grc_1.write('\n')

	for entry in haplotype_data:
		sample_ids = list(entry.keys())
		for id in sample_ids:
			data = entry[id]
			grc_1.write(id+'\t')
			for gene in drug_resistance_loci['genes']:
				grc_1.write(data[gene] + '\t')
			grc_1.write('\n')

if __name__ == "__main__":
	main()
