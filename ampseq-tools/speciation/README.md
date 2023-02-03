# Speciation Module
A module to perform speciation based on the current amplicon production pipeline.

## Workflow - per sample
![speciation_flowchart](./speciation_flow.png)

## Requirements
* Python > 3.8
* pandas
* tqdm

# Run speciate on batch
```
<ampseq_pipeline_repo_location>/ampseq-tools/speciation/speciate \
"<path_to_input_genotype_files>" \
<path_to_barcode_file> \
<path_to_config_file>
```
# Run unit tests
```
<ampseq_pipeline_repo_location>/ampseq-tools/speciation/test
```