from AminoAcidCalling import AminoAcidCaller

def main():

	#rudimentary implementation to run methods only
	



	a = AminoAcidCaller("SAMPLE_ID", "./genotype_file.tsv", './config.json')
	a.call_haplotypes()


if __name__ == "__main__":
	main()

	
