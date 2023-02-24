# Barcoding Module

## Help Message
```
>> ./barcoding --help

usage: Amplicon Barcode Production [-h] [--outfile OUTFILE] [--pbar]
                                   [--ncpus NCPUS] [--sample_col SAMPLE_COL]
                                   genotype_files config

A package to perform barcode production based on the production amplicon
pipeline

positional arguments:
  genotype_files        Path to input genotype file(s)
  config                Path to config json file

optional arguments:
  -h, --help            show this help message and exit
  --outfile OUTFILE     Path to directory to output results (default:
                        barcode_results.txt)
  --pbar                Show a progress bar while running
  --ncpus NCPUS         No. cpus to use in processing
  --sample_col SAMPLE_COL
                        Column to access for sample ID information in genotype
                        file (default: SampleID)
```

## Requirements
* Python>3.7
* tqdm