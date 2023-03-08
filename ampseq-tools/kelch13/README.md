# Kelch13 Mutation Calling Module

## Help Message
```
>> ./grc_kelch13_mutation_caller.py --help

usage: grc_kelch13_mutation_caller.py [-h] [--genotype_files GENOTYPE_FILES]
                        [--config CONFIG] [--output_file OUTPUT_FILE]
                        [--kelch_reference_file KELCH_REFERENCE_FILE]
                        [--codon_key_file CODON_KEY_FILE]

A package to perform kelch13 mutation calling based on the production amplicon
pipeline

arguments:
  -h, --help            show this help message and exit
  --genotype_files GENOTYPE_FILES
                        Path to input genotype file(s)
  --config CONFIG       Path to config json file
  --output_file OUTPUT_FILE
                        Path to directory to output results (default:
                        kelch13_mutation_calls.txt)
  --kelch_reference_file KELCH_REFERENCE_FILE
                        Path to kelch13 reference file
  --codon_key_file CODON_KEY_FILE
                        Path to codon key file
```

# Run kelch13 mutation calling on batch
```
<ampseq_pipeline_repo_location>/ampseq-tools/kelch13/grc_kelch13_mutation_caller.py \
--genotype_files "<path_to_input_genotype_files>" \
--config <path_to_config_file> \
--kelch_reference_file <path_to_kelch_reference_file> \
--codon_key_file <path_to_codon_key_file>
```

## Requirements
* Python>3.7
