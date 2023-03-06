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
usage: grc_speciate [-h] [--genotype_files GENOTYPE_FILES]
                    [--barcodes_file BARCODES_FILE] [--config CONFIG]
                    [--output_file OUTPUT_FILE] [--pbar] [--ncpus NCPUS]
                    [--output_debug_path OUTPUT_DEBUG_PATH]

A package to perform speciation based on the production amplicon pipeline

optional arguments:
  -h, --help            show this help message and exit
  --genotype_files GENOTYPE_FILES
                        Input list of all genotype files (TSV format) to be
                        run through the speciation program
  --barcodes_file BARCODES_FILE
                        Path to barcodes output file for querying
  --config CONFIG       Path to config json
  --output_file OUTPUT_FILE
                        Path to output file
  --pbar                Show a progress bar while running
  --ncpus NCPUS         No. cpus to use in processing
  --output_debug_path OUTPUT_DEBUG_PATH
                        Ouput directory to save debug files. If not provided
                        files will not be output.
```