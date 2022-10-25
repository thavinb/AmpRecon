# AmpSeq-pipeline

Ampseq is a bioinformatics analysis pipeline for amplicon sequencing data. Currently supporting alignment and SNP variant calling on paired-end Illumina sequencing data.

The pipeline has capabilities to generate variant call format (VCF) files directly from a bcl directory, as well as starting from aligned CRAM formatted files stored in Sanger's internal file storage system (iRODS). Our pipeline allows configurable reference handling, allowing high-throughput processing of data against multiple amplicon panels in a single pipeline run.

We opted for industry standard Burrows-Wheeler Aligner (BWA) coupled to GATK4's haplotypecaller and genotypegvcf to call variants.  

# Pipeline summary
Using the default run options, ampseq performs the following tasks:
- Converts a .bcl directory into a BAM formatted file ([bambi i2b](https://wtsi-npg.github.io/bambi/#i2b))
- Decodes multiplexed BAM file ([bambi decode](https://wtsi-npg.github.io/bambi/#decode))
- Sequencing adapter contamination removal ([biobambam2 bamadapterfind](https://manpages.ubuntu.com/manpages/focal/en/man1/bamadapterfind.1.html))
- BAM to CRAM conversion ([samtools split](http://www.htslib.org/doc/samtools-split.html))
- Data cleaning and QC processes prior to alignment and SNP calling ([cram to bam](https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline/-/blob/develop/workflows/pipeline-subworkflows/cram-to-bam.nf))
- Alignment to full reference / amplicon panel ([realignment](https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline/-/blob/develop/workflows/pipeline-subworkflows/realignment.nf))
- SNP calling using GATK4 tools haplotypecaller and genotypegvcf joint genotyping protocol ([genotyping](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4))

# Quick-start guide to Ampseq

Ampseq v0.0.1 is currently only available on Sanger's internal high performance compute (HPC) clusters with access to nfs storage.

**1. Setup**

clone the repository
```
git clone https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline.git
cd ./ampseq-pipeine/
``` 

**2. Run**

Load nextflow module
```
module load nextflow/22.04.0-5697
```

To run from the **in-country** entry point:

```
nextflow ../path/to/ampseq-pipeline/main.nf -profile sanger_lsf \
                --execution_mode in-country --run_id 21045 \
                --bcl_dir /path/to/my_bcl_dir/ --lane 1 \
                --study_name test --read_group rg_test --library lib
```

To run from the **iRODS** entry point:

```
nextflow ../ampseq-pipeline/main.nf -profile sanger_lsf \ 
        --execution_mode irods \ 
        --irods_manifest ./input/irods_smallset.tsv
```

Use `-profile sanger_lsf` to make nextflow be able to submit task to the farm lsf queue.
Use `-profile sanger_default` to run on the farm but local (this should be used only for dev porpuse).

## Parameters

Absolutely required
```
exectution_mode : sets the entry point for the pipeline ("irods" or "in-country")
```

Required for **in-country**
```
run_id : id to be used for the batch of data to be processed
bcl_dir: path to a miseq directory
```
The pipeline needs values for `--lane , --study_name, --read_group, --library` are requested, but the content is irrelevant for the pipeline execution.
Therefore, choose wisely.

Required for **iRODS**
```
irods_manifest : an tsv containing information of irods data to fetch
```
An example of an irods manifest tsv is provided at [add path to example]

To use **S3**
```
s3_uuid : <str> a universally unique id which will be used to fetch data from s3,
s3_bucket_input: <str> s3 bucket name to fetch data from, if is not provided, the pipeline will not retrieve miseq runs from s3

upload_to_s3: <bool> sets if needs to upload output data to an s3 bucket
s3_bucket_output: <str> s3 bucket name to upload data to
```
### iRODS Manifest
The iRODS must be a `.tsv` and the pipeline expects to find the following columns headers:

* `sample_id`: a sample identification "tag";

* `primer_panel`: primer pannel name to be used (must match exactly what is provided at `pannel_name` of the `pannels_settings.csv`);

* `irods_path`: full valid iRODS path for a `.cram` file (ex: `/seq/illumina/runs/38/38344/lane2/plex1/38344_2#1.cram`).

The `.tsv` may have more columns at any order, but those are the only ones which will be considered.
The pipeline builds an "internal id" set as `<cram_filename>_<sample_id>_<primer_panel>`, therefore, the pipeline will check if any combination of those values at the manifest are unique. If not, an error will be raised and the pipeline run will stop.
An example of a valid manifest can be found at this repository (`test_data/irods_mnf.tsv`).

### Pannel Settings
The ampseq pipeline relies on a `pannels_settings.csv` to define which files it should use on key pipeline steps according to the pannel name provided for a given sample.
Currently, this `.csv` should look like the example below:

```
pannel_name,aligns_to,maps_to_regions_of
PFA_GRC1_v1.0,/path/to/pannels_resources/grc1/,/path/to/PFA_GRC1_v1.0.annotation.regions.txt
PFA_GRC2_v1.0,/path/to/pannels_resources/grc2/,/path/to/PFA_GRC2_v1.0.annotation.regions.txt
PFA_Spec,/path/to/pannels_resources/spec/,/path/to/PFA_Spec_v1.0.annotation.regions.txt
```

* `pannel_name` : Defines the string it should look for a given pannel, this strings should be the same provided by the user (via samplesheet or irods_manifest).

* `aligns_to` : Defines which directory it should look to get the `.fasta` (and associated index files) to use for the alignment, namely  `IN_COUNTRY:CRAM_TO_BAM:align_bam` and `COMMON:REALIGNMENT:aligns_bam`.

* `maps_to_regions_of` : Defines which annotation file should use for the `COMMON:REALIGNMENT:read_count_per_region`.

This pannel settings system aims to detach the experimental design from the inner works of the pipeline and make it easier to experiment with its key steps. A custom `.csv` can be set to the pipeline by using the flag `--pannels_settings`. If the user does not provide a `--pannels_settings`, the pipeline default behaviour is to rely on files available at the repo (check `pannels_resources` dir).

### Genotyping Settings
The genotyping portion of the ampseq pipeline uses a variety of parameters defined in a `methods.config` file.
The following parameters should be present within this file:

```
gatk3: <str> path to GATK3 GenomeAnalysisTK.jar file.
combined_vcf_file1: <str> known SNPs database file. Used to prevent BaseRecalibrator from using regions surrounding polymorphisms.
combined_vcf_file2: <str> known SNPs database file. Used to prevent BaseRecalibrator from using regions surrounding polymorphisms.
combined_vcf_file3: <str> known SNPs database file. Used to prevent BaseRecalibrator from using regions surrounding polymorphisms.
conserved_bed_file: <str> file containing genomic intervals the GATK BaseRecalibrator command operates over in the bqsr.nf process.
gatk_base_recalibrator_options: <str> input settings containing the supplied known sites files paths and intervals file path for the BaseRecalibrator command in the bqsr.nf process.
alleles_fn: <str> file containing genomic intervals the GATK GenotypeGVCFs command operates over in the genotype_vcf_at_given_alleles.nf process.
```

### Containers

By default, the pipeline will look for the containers at `/nfs/gsu/team335/ampseq-containers/`. Another directory to look for the containers can be set using the `--containers_dir` flag at the nextflow command line.
All recipes for the ampseq containers can be found at the `containers/` directory of this repository. Use the following command to build it:

```
cd /path/to/ampseq-pipeline/containers/
bash buildContainers.sh
```

The building process take a few minutes to finish and all necessary `.sif` files to run the pipeline will be generated.

## Support
[who should someone talk to regarding the maintenance and usage of the pipeline]

## Authors and acknowledgment
Show our appreciation to those who have contributed to the project.

## License
[add a license]

## Project status
[prototype]
