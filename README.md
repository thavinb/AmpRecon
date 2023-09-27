# AmpSeq

Ampseq is a bioinformatics analysis pipeline for amplicon sequencing data. Currently supporting alignment and SNP variant calling on paired-end Illumina sequencing data.

The pipeline has capabilities to output [Genetic Report Cards (GRCs)](https://www.malariagen.net/sites/default/files/GRC_UserGuide_10JAN19.pdf) and readcounts per pannel files  directly from [Binary Base Calls (BCLs)](https://emea.illumina.com/informatics/sequencing-data-analysis/sequence-file-formats.html) files, as well as starting from aligned [CRAM](https://www.sanger.ac.uk/tool/cram/) formatted files stored in Sanger's internal file storage system which is based on [iRODS](https://irods.org) (Integrated Rule-Oriented Data System). In addition, the pipeline also outputs [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) per lanelet and [VCFs](https://samtools.github.io/hts-specs/VCFv4.2.pdf) per sample. Our pipeline allows configurable reference handling, allowing high-throughput processing of data against multiple amplicon panels in a single pipeline run. We allow both plasmodium falciparum and vivax species to be run through the invocation of a configuration file in the conf sub-directory of this repository


# Quick-start guide to Ampseq

Ampseq v0.0.1 is currently only available on Sanger's internal high performance compute (HPC) clusters with access to nfs storage.

## 1. Setup

### 1.1 - requirements

The pipeline requires:

-  [Nextflow](https://github.com/nextflow-io/nextflow), the pipeline was developed and tested on [version 22.04](https://github.com/nextflow-io/nextflow/releases/tag/v22.04.4). 

- [Singularity](https://github.com/sylabs/singularity), the pipeline was developed and tested on the singularity [version 3.6.4](https://github.com/apptainer/singularity/releases/tag/v3.6.4).

> **WARNING**
> If Singularity is not available, the pipeline assumes all the softwares with the correct versions are available on the execution environment.

### 1.2 - Download code base and build containers

All recipes for the ampseq containers can be found at the `containers/` directory of this repository, The building process take a few minutes to finish and all necessary `.sif` files to run the pipeline will be generated on the same dir.

```
# clone repo
git clone --recurse-submodules https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline.git
# build containers
cd ./ampseq-pipeine/containers/
bash buildContainers.sh
```

---
## 2 - Run the pipeline 

Assuming nextflow is available at the command line, to run from the **in-country** entry point:

```
nextflow /path/to/ampseq-pipeline/main.nf -profile sanger_lsf \
                --execution_mode in-country \
                --run_id 21045 \
                --bcl_dir /path/to/my_bcl_dir/ --ena_study_name test --manifest_path manifest.tsv \
                --containers_dir /path/to/containers_dir/
                -c /path/to/species/config
```

To run from the **iRODS** entry point:

```
nextflow /path/to/ampseq-pipeline/main.nf -profile sanger_lsf \
        --execution_mode irods \
        --run_id 21045 \
        --irods_manifest /path/to/irods_manifest.tsv
        --containers_dir ./containers_dir/
        -c /path/to/species/config 
```


To run from the **FASTQ** entry point:

```
nextflow /path/to/ampseq-pipeline/main.nf -profile sanger_lsf \
        --execution_mode fastq \
        --run_id 21045 \
        --fastq_manifest /path/to/fastq_manifest.tsv
        --containers_dir ./containers_dir/
        -c /path/to/species/config 
```


> **NOTE**
> If **on the farm**, Nextflow can be made available by loading its module.
>```
>module load nextflow/22.04.0-5697
>```

Use `-profile sanger_lsf` to make nextflow be able to submit task to the farm lsf queue.
Use `-profile sanger_default` to run on the farm but local (this should be used only for dev porpuse).
To use a panel resources different than the ones provided at this repository, the user needs to provide a custom pannels settings csv via `--panels_settings`.

---

## Parameters

Absolutely required

```
execution_mode : sets the entry point for the pipeline ("irods" or "in-country")
species_config/-c : stages the relevant reference files to run the species selected. Two configuration files for plasmodium vivax and plasmodium falciparum are provided in the repository in the conf sub-directory
```

Required for **in-country**

```
run_id : id to be used for the batch of data to be processed. Also added as a prefix to the output GRC files
bcl_dir: path to a miseq directory
manifest_path: path to the manifest tab separated values file.
```

The pipeline needs a value for `--ena_study_name` to be supplied, but the content is irrelevant for the pipeline execution.
Therefore, choose wisely.

Required for **iRODS**

```
irods_manifest : an tsv containing information of irods data to fetch
run_id : Added as a prefix to the output GRC files
```

An example of an irods manifest tsv is provided at [add path to example]

Required for **FASTQ**

```
fastq_manifest : an tsv containing information for the unpaired fastq data
run_id : Added as a prefix to the output GRC files
```

An example of an fastq manifest tsv is provided at [add path to example]. Currently only unpaired fastqs are supported.

To use **S3**

```
upload_to_s3: <bool> sets if needs to upload output data to an s3 bucket
s3_bucket_output: <str> s3 bucket name to upload data to
```


**Genotyping**

The genotyping portion of the ampseq pipeline requires a SNPs annotation VCF file for each reference genome used. This file location must be set at the `snp_list` column of the `panel_settings.csv`.

**Containers**

By default, the pipeline will look for the containers at `/nfs/gsu/team335/ampseq-containers/`. A different directory to look for the containers can be set using the `--containers_dir` flag at the nextflow command line.

---
## Input Files
### In Country Manifest

The in country manifest file must be a `.tsv` and the pipeline expects to find the following columns headers:

- `sample_id`: a sample identifier;

- `primer_panel`: primer panel name to be used (must match exactly what is provided at `panel_name` of the `panels_settings.csv`);

- `barcode_number`: a unique number for each lanelet;

- `barcode_sequence`: two DNA barcode sequences separated by a hyphen;

- `partner_sample_id`: name allocated to the sample by partner. Metadata to be added to the final GRC files

- `collection_date`: sample collection date. Metadata to be added to the final GRC files

- `collection_location`: name of the location within the country where the sample was collected. Metadata to be added to the final GRC files

- `collection_country`: name of country the sample was collected in. Metadata to be added to the final GRC files

- `study`: full study ID of the sample. Metadata to be added to the final GRC files

- `well`: a well identifier

- `plate_name`: a plate identifier

```
sample_id	primer_panel	barcode_number	barcode_sequence	partner_sample_id	collection_date	collection_location	collection_country	study	well	plate_name
ILCM4453	PFA_GRC1_v1.0	1	ATCACGTT-GTACTGAC	ILCM4453	2021-07-16	05at3a Samroung Romdul health center	Cambodia	130as-ICS	A01	PLATE_RCN_00190
RAD3CC01	PFA_GRC2_v1.0	2	CGATGCAT-GTACTACC	TTCN3A01	2021-09-12	05a0qt Chambak health center	Cambodia	130as-ICS	A02	PLATE_RCN_00190
RIALMC99	PFA_Spec	3	TTAACACT-GTACTGAC	IL21939L	2021-10-21	0406xq Chambak health center	Cambodia	130as-ICS	A03	PLATE_RCN_00190
```

### iRODS Manifest

The iRODS manifest file must be a `.tsv` and the pipeline expects to find the following columns headers:

- `sample_id`: a sample identification "tag", which is used on the pipeline output file names;

- `primer_panel`: primer panel name to be used (must match exactly what is provided at `panel_name` of the `panels_settings.csv`);

- `irods_path`: full valid iRODS path for a `.cram` file (ex: `/seq/illumina/runs/38/38344/lane2/plex1/38344_2#1.cram`).

- `partner_sample_id`: name allocated to the sample by partner. Metadata to be added to the final GRC files

- `collection_date`: sample collection date. Metadata to be added to the final GRC files

- `collection_location`: name of the location within the country where the sample was collected. Metadata to be added to the final GRC files

- `collection_country`: name of country the sample was collected in. Metadata to be added to the final GRC files

- `study`: full study ID of the sample. Metadata to be added to the final GRC files

The `.tsv` may have more columns at any order, but those are the only ones which will be considered.
The pipeline builds an "internal id" set as `<cram_filename>_<sample_id>_<primer_panel>`, therefore, the pipeline will check if any combination of those values at the manifest are unique. If not, an error will be raised and the pipeline run will stop.
The content of valid manifest should look like the example bellow:

```
irods_path	sample_id	primer_panel	study_name	pipeline_id	taxon_id	common_name	name	supplier_name	donor_id	instrument_model	qc_complete	id_run	lane	tag	qc	WG_lane
/seq/29632/29632_1#55.cram	ILL411270	PFA_GRC1_v1.0	Team 112 R&D	GBS	5833	Plasmodium Falciparum	3429STDY7977888	RCN15139	3429STDY7977888	MiSeq	2019-05-30 03:38:57	29632	1	55	1	29632_1#55
/seq/29632/29632_1#149.cram	LMLPP1571	PFA_GRC2_v1.0	Team 112 R&D	GBS	5833	Plasmodium Falciparum	3429STDY7977888	RCN15139	3429STDY7977888	MiSeq	2019-05-30 03:38:57	29632	1	149	1	29632_1#149
/seq/26381/26381_1#808.cram	JHG3639016I	PFA_Spec	Team 112 R&D	GBS	5833	Plasmodium Falciparum	3429STDY7977859	RCN15110        3429STDY7977888	MiSeq	2019-05-30 03:38:57	29632	1	149	1	29632_1#256
```
### FASTQ Manifest

The FASTQ manifest file must be a `.tsv` and the pipeline expects to find the following columns headers:

- `sample_id`: a sample identification "tag", which is used on the pipeline output file names;

- `primer_panel`: primer panel name to be used (must match exactly what is provided at `panel_name` of the `panels_settings.csv`);

- `fastq_path`: full valid FASTQ path for a `.fastq` file.

```
sample_id	primer_panel	fastq_path
sample1	PFA_GRC1_v1.0	/path/to/fastq/sample1.fastq
sample2	PFA_GRC2_v1.0	/path/to/fastq/sample2.fastq
sample3	PFA_Spec	/path/to/fastq/sample3.fastq
```


### Species configuration file

Ampseq relies on a configuration file that specifies the species specific run settings. To invoke the pipeline it is imperative that one of these configuration files is passed at the command line using nextflow run's -c flag. If the user does not have a configuration file. They may use one of the ones provided in the conf sub-directory. The species configuration file points to files present at the [ampliconresources submodule](https://gitlab.internal.sanger.ac.uk/malariagen1/ampliconresources/-/tree/main/).

### Panel Settings

The ampseq pipeline relies on a `panels_settings.csv` to define which files it should use on key pipeline steps according to the panel name provided for a given sample.
Currently, this `.csv` should look like the example below:

```
panel_name,reference_file,design_file,snp_list
PFA_GRC1_v1.0,/path/to/PFA_GRC1_v1.0.fasta,/path/to/PFA_GRC1_v1.0.regions.txt,/path/to/PFA_GRC1_v1.0.annotation.vcf
PFA_GRC2_v1.0,/path/to/PFA_GRC2_v1.0.fasta,/path/to/PFA_GRC2_v1.0.regions.txt,/path/to/PFA_GRC2_v1.0.annotation.vcf
PFA_Spec,/path/to/PFA_Spec.fasta,/path/to/PFA_Spec.regions.txt,/path/to/PFA_Spec.annotation.vcf
```

- `panel_name` : Defines the string it should look for a given panel, this strings should be the same provided by the user (via the supplied manifest file).

- `reference_file` : Path to the `reference.fasta` for use in alignment. Reference index files (.fai, .amb, .ann, .bwt, .pac and .sa) and a sequence dictionary file (reference_file_name.dict) should also be found at this location.

- `design_file` : Defines which annotation file should use for the `COMMON:REALIGNMENT:read_count_per_region`.

- `snp_list` : Path to the SNP list file, used as a targets file for BCFtools mpileup.

The aim of this panel settings system is to detach the experimental design from the inner works of the pipeline and make it easier to experiment with its key steps. A `.csv` **must be privded** to the pipeline via `--panels_settings`.

### GRC Creation Requirements

The GRC creation of the pipeline requires several additional files. These include a GRC settings JSON file (`--grc_settings_file_path`), a chromKey file (`-- --chrom_key_file_path`), a codonKey file (`--condon_key_file_path`), a Kelch13 reference sequence file (`--kelch_reference_file_path`) and a drug resistance loci information file (`--drl_information_file_path`), which must be provided via nextflow parameters.

---

## Running NF-tests

The unit tests for the workflow are implemented using [NF-test](https://code.askimed.com/nf-test/).
If not available already on the CLI,
### Install NF-test

1. Download NF-test:

```{bash}
wget -qO- https://code.askimed.com/install/nf-test | bash
```

2. Create an alias for the file you just downloaded, be sure you can execute the file.

```{bash}
alias nf-test="/path/to/my/nf-test"
```

### Run tests

On the repository directory, run:

```{bash}
nf-test test tests/workflows/sanger_irods_to_reads.nf.test --profile sanger_default
nf-test test tests/workflows/miseq_to_reads.nf.test --profile sanger_default
```
