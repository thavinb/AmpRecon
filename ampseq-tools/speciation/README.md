# Speciation Module
A module to perform speciation based on the current amplicon production pipeline.

## Workflow - per sample
![speciation_flowchart](./speciation_flow.png)

## Requirements
* Python > 3.8
* tqdm

# Run speciate on batch
```
<ampseq_pipeline_repo_location>/ampseq-tools/speciation/grc_speciation.py \
--genotype_files "<path_to_input_genotype_files>" \
--barcodes_file <path_to_barcode_file> \
--config <path_to_config_file>
```

# Run unit tests
```
<ampseq_pipeline_repo_location>/ampseq-tools/speciation/test
```

# Help Message
```
usage: speciate [-h] [--outfile OUTFILE] [--pbar]
                [--ncpus NCPUS]
                [--output_debug_path OUTPUT_DEBUG_PATH]
                [--sample_col SAMPLE_COL]
                [--chrom_regex CHROM_REGEX]
                input_genotype_files barcodes config

A package to perform speciation based on the production amplicon pipeline

positional arguments:
  input_genotype_files  Input list of all genotype files (TSV format) to be
                        run through the speciation program
  barcodes              Path to barcodes output file for querying
  config                Path to config json (default: config.json)

optional arguments:
  -h, --help            show this help message and exit
  --outfile OUTFILE     Path to output file (default: ./barcodes.tsv)
  --pbar                Show a progress bar while running
  --ncpus NCPUS         No. cpus to use in processing
  --output_debug_path OUTPUT_DEBUG_PATH
                        Ouput directory to save debug files. If not provided
                        files will not be output.
  --sample_col SAMPLE_COL
                        Column to access for sample ID information in genotype
                        file (default: SampleID)
  --chrom_regex CHROM_REGEX
                        A regex for matching on spec panel chromosomes
                        (default: "^Spec_[12]_(falciparum|vivax)$")
```