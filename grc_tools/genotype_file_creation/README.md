# Genotype File Writing Module

## Help Message
```
>> ./write_genotypes_file.py --help

usage: write_genotypes_file.py [-h] [--vcf_files VCF_FILES]
                        [--output_file OUTPUT_FILE] [--sample_id SAMPLE_ID]
                        [--chromKey_file CHROMKEY_FILE_PATH]

A package to perform genotype file writing based on the production amplicon
pipeline

optional arguments:
  -h, --help            show this help message and exit
  --chromosome_column_name CHROMOSOME_COLUMN_NAME
                        Name of the column in the chromKey file to try to match chromosome to
  --locus_column_name LOCUS_COLUMN_NAME
                        Name of the column in the chromKey file to try to match position to
  --min_total_depth MIN_TOTAL_DEPTH
                        The number of reads at a record must exceed this value
  --het_min_allele_depth HET_MIN_ALLELE_DEPTH
                        The number of reads for a particular allele at a heterozygous record must exceed this value
  --het_min_allele_proportion HET_MIN_ALLELE_PROPORTION
                        The proportion of reads for particular allele at a heterozygous record must exceed this value
```

# Run genotype file writing script on batch
```
<ampseq_pipeline_repo_location>/ampseq-tools/genotype_file_creation/write_genotypes_file.py \
--vcf_files <vcf_file_path_list> \
--output_file <name_of_output_file> \
--sample_id <sample_identifier> \
--chromKey_file <path_to_chromKey_file> \
--chromosome_column_name <chromKey_chromosome_column_name> \
--locus_column_name <chromKey_locus_column_name> \
--min_total_depth <minimum_depth_value_for_locus> \
--het_min_allele_depth <minimum_depth_value_for_allele> \
--het_min_allele_proportion <minimum_proportion_value_for_allele>
```

## Requirements
* Python>3.7
* PyVCF>0.6.8
