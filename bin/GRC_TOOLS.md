# Documentation for Python Scripts Used in GRC Creation  

**NB**: All scripts described here require Python>=3.8; other requirements are noted as applicable.

- [Genotype File Creation](#genotype-file-creation)  
- [Kelch13](#kelch13-mutation-detection)  
- [Plasmepsin CNV Detection](#plasmepsin-copy-number-variation-detection)  
- [Barcode Generation](#barcode-generation)  
- [Species Detection](#species-detection)  
- [Complexity-of-Infection Estimation](#complexity-of-infection-estimation)  

# Genotype File Creation  

**Script** : `write_genotypes_file.py`  
**Description**: Generate genotype file as part of GRC generation  
**Called From**: `grc_assemble_genotype_file.nf`  

#### Requirements

* PyVCF>0.6.8  

## Usage  

```
write_genotypes_file.py [-h] [--vcf_files VCF_FILES]
                        [--output_file OUTPUT_FILE] [--sample_id SAMPLE_ID]
                        [--chromKey_file CHROMKEY_FILE_PATH]
```  

## Parameters  
### Required  

- `-h, --help` : Print help and exit.  
- `--vcf_files` (paths): List of VCF files to use as input for genotype file generation.  
- `--output_file` (path) : Path to output genotype file.  
- `--sample_id` (str) : Sample identifier.  
- `--chromKey_file` (path) : Path to chromKey file.  

### Optional
- `--chromosome_column_name` : Name of the column in the chromKey file to try to match chromosome to.  
- `--locus_column_name` : Name of the column in the chromKey file to try to match position to.  
- `--min_total_depth` : Minimum depth of coverage required to consider a record.  
- `--het_min_allele_depth` : Minimum depth of coverage required to consider a _heterozygous_ record.  
- `--het_min_allele_proportion` : Minimum heterozygous allele proportion required to consider a record.  

[**(&uarr;)**](#documentaion-for-python-scripts-used-in-grc-creation)  

---

# Kelch13 Mutation Detection  

**Script** : `grc_kelch13_mutation_caller.py`  
**Description**: Call kelch13 mutations from input genotype file.  
**Called From**: `grc_kelch13_mutation_caller.nf`  

## Usage  

```
grc_kelch13_mutation_caller.py [-h] [--genotype_files GENOTYPE_FILES]
                        [--config CONFIG] [--output_file OUTPUT_FILE]
                        [--kelch_reference_file KELCH_REFERENCE_FILE]
                        [--codon_key_file CODON_KEY_FILE]
```  

## Parameters  
### Required  

- `-h, --help` : Print help and exit.  
- `--genotype_files` (path) : Path to input genotype file(s).  
- `--config` (path) : Path to config json file.  
- `--output_file` (path) [Default: kelch13_mutation_calls.txt] : Path to directory to output results  
- `--kelch_reference_file` (path) : Path to kelch13 reference file.  
- `--codon_key_file` (path) : Path to codon key file.  

[**(&uarr;)**](#documentaion-for-python-scripts-used-in-grc-creation)  

---

# Plasmepsin Copy-number Variation Detection  

**Script** : `grc_plasmepsin_cnv_caller.py`  
**Description**: Call copy-number variants for the plasmepsin locus.  
**Called From**: `grc_plasmepsin_cnv_caller.nf`  

## Usage  

```
grc_plasmepsin_cnv_caller.py [-h] [--genotype_files GENOTYPE_FILES]
                        [--config CONFIG] [--output_file OUTPUT_FILE]
```  

## Parameters  
### Required  

- `-h, --help` : Print help and exit.  
- `--genotype_files` (path): Path to input genotype file(s)  
- `--config CONFIG` (path): Path to config json file  
- `--output_file` (path) [Default: plasmepsin_variant_calls.txt] : Path to directory to output results  

[**(&uarr;)**](#documentaion-for-python-scripts-used-in-grc-creation)  

---

# Barcode Generation  

**Script** : `grc_barcoding.py`  
**Description**: Generate barcodes for GRC.  
**Called From**: `grc_barcodeing.nf`

#### Requirements  

* tqdm  

## Usage  

```
grc_barcoding.py [-h] [--genotype_files GENOTYPE_FILES]
                        [--config CONFIG] [--output_file OUTPUT_FILE] [--pbar]
                        [--ncpus NCPUS]
```

## Parameters  
### Required  

- `-h, --help` : Print help and exit.  
- `--genotype_files` (path): Path to input genotype file(s)  
- `--config CONFIG` (path): Path to config json file  
- `--output_file` (path) [Default: barcode_results.txt] : Path to directory to output results  

### Optional  

- `--pbar` : Show a progress bar while running  
- `--ncpus` (int): Number of SPUs to use in  

[**(&uarr;)**](#documentaion-for-python-scripts-used-in-grc-creation)  

---

# Species Detection  

**Script** : `grc_speciate.py`  
**Description**: Run species detection on sample data  
**Called From**: `grc_speciate.nf`  

#### Requirements  

* tqdm  

## Usage  

```
grc_speciate.py [-h] [--genotype_files GENOTYPE_FILES]
                    [--barcodes_file BARCODES_FILE] [--config CONFIG]
                    [--output_file OUTPUT_FILE] [--pbar] [--ncpus NCPUS]
                    [--output_debug_path OUTPUT_DEBUG_PATH]
```  

## Parameters  
### Required  

- `-h, --help` : Print help and exit.  
- `--genotype_files` (paths): List of all genotype files (TSV format) to be run through the speciation program  
- `--barcodes_file` (path): Path to barcodes output file for querying  
- `--config` (path): Path to config json  
- `--output_file` (path): Path to output file  

### Optional  

- `--pbar` : Show a progress bar while running  
- `--ncpus` (int): Number of SPUs to use in  
- `--output_debug_path` (path): Ouput directory to save debug files. If not provided files will not be output.

[**(&uarr;)**](#documentaion-for-python-scripts-used-in-grc-creation)  

---

# Complexity-of-Infection Estimation  

**Script** : `grc_process_mccoil_io.py`  
**Description**: Generate input files for, and subsequently run, THEREALMcCOIL
**Called From**: `grc_estimate_coi.nf`

#### Requirements  

* [THEREALMcCOIL - (slightly modified version)](https://github.com/AMarinhoSN/THEREALMcCOIL)
* Singularty

## Usage  

```
grc_process_mccoil_io.py [-h] [-write_mccoil_in] [-write_coi_grc] [--barcodes_files BARCODES_FILE(s)]
                            [--barcode_def_file BARCODE_DEF_FILE]
                            [--mccoil_sum_file MCCOIL_SUM_FILE(s)]
                            [--output_file OUTPUT_FILE]
```  

Ideally, use Singularity to build from the definition file in THEREALMcCOIL repository - this will handle all dependencies. Within the container, typical usage would be as follows (see the repository for more information on the options taken by THEREALMcCOIL):  

```
## generate McCOIL input files
python3 grc_process_mccoil_io.py -write_mccoil_in \
            --barcodes_files $BARCODES_IN --config $BARCODE_DEF \
            --output_file RUN_ID.tsv

## run McCOIL
Rscript /app/THEREALMcCOIL/runMcCOIL.R -i RUN_ID.tsv \
            --totalRun NTOTAL --totalBurnIn NBURN --seed 123456 \
            --outPrefix RUN_ID --maxCOI 20 --M0 5

## convert McCOIL outputs to GRC-ready form
python3 grc_process_mccoil_io.py -write_coi_grc \
            --mccoil_sum_file RUN_ID_summary.txt \
            --output_file RUN_ID_out.grc
```  

## Parameters  
### Required  

- `-h, --help` : Print help and exit.  
- `-write_mccoil_in` : Prepare McCOIL input files from barcodes  
- `-write_coi_grc` : Prepare McCOIL input files from barcodes  
- `--barcodes_files` (paths): Path(s) to barcode `.tsv` file for a batch of samples  
- `--barcode_def_file` (path): Path to a json file with barcodes definitions  
- `--mccoil_sum_file` (path): Path to McCOIL summary output file  
- `--output_file` (path) [Default:  <./McCOIL_in.tsv | ./coi.grc >]: Path for the McCOIl input or GRC file to be written. If more than one input file is provided, the corresponding outputs will have the same basename as the respective input files prefixed to the standard file suffixes.  

### Note

Inputs and outputs of the original The REAL McCOIL should work just fine, however, here we use a [custom version](https://github.com/AMarinhoSN/THEREALMcCOIL). Small tweaks were necessary to containarize this software and improve usability by adding a command-line interface.  

[**(top &uarr;)**](#documentaion-for-python-scripts-used-in-grc-creation)  
