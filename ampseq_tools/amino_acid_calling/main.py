from AminoAcidCalling import AminoAcidCaller

def main():

	#rudimentary implementation to run methods only
	
	a = AminoAcidCaller("test_sample", "./genotype_file.tsv")
	a.call_haplotypes()



if __name__ == "__main__":
	main()
