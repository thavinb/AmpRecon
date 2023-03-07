# Barcoding Module

## Help Message
```
>> ./grc_barcoding.py --help

usage: grc_barcoding.py [-h] [--genotype_files GENOTYPE_FILES]
                        [--config CONFIG] [--output_file OUTPUT_FILE] [--pbar]
                        [--ncpus NCPUS]

A package to perform barcode production based on the production amplicon
pipeline

optional arguments:
  -h, --help            show this help message and exit
  --genotype_files GENOTYPE_FILES
                        Path to input genotype file(s)
  --config CONFIG       Path to config json file
  --output_file OUTPUT_FILE
                        Path to directory to output results (default:
                        barcode_results.txt)
  --pbar                Show a progress bar while running
  --ncpus NCPUS         No. cpus to use in processing
```

# Run speciate on batch
```
<ampseq_pipeline_repo_location>/ampseq-tools/barcode/grc_barcoding.py \
--genotype_files "<path_to_input_genotype_files>" \
--config <path_to_config_file>
```

## Requirements
* Python>3.7
* tqdm