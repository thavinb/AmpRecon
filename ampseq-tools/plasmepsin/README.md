# Plasmepsin CNV Calling Module

## Help Message
```
>> ./grc_plasmepsin_cnv_caller.py --help

usage: grc_plasmepsin_cnv_caller.py [-h] [--genotype_files GENOTYPE_FILES]
                        [--config CONFIG] [--output_file OUTPUT_FILE]

A package to perform plasmepsin copy number variation calling based on the production amplicon
pipeline

arguments:
  -h, --help            show this help message and exit
  --genotype_files GENOTYPE_FILES
                        Path to input genotype file(s)
  --config CONFIG       Path to config json file
  --output_file OUTPUT_FILE
                        Path to directory to output results (default:
                        plasmepsin_variant_calls.txt)
```

# Run plasmepsin copy number variation calling on batch
```
<ampseq_pipeline_repo_location>/ampseq-tools/plasmepsin/grc_plasmepsin_cnv_caller.py \
--genotype_files "<path_to_input_genotype_files>" \
--config <path_to_config_file>
```

## Requirements
* Python>3.7
