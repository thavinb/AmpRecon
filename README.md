# AmpSeq-pipeline

Ampseq is a bioinformatics analysis pipeline for amplicon sequencing data. Currently supporting alignment and SNP variant calling on paired-end Illumina sequencing data.

The pipeline has capabilities to generate variant call format (VCF) files directly from a bcl directory, as well as starting from aligned CRAM formatted files stored in Sanger's internal file storage system (iRODS).  

# Pipeline summary
Using the default run options, ampseq performs the following tasks:
- Converts a .bcl directory into a BAM formatted file ([bambi i2b](https://wtsi-npg.github.io/bambi/#i2b))
- Decodes multiplexed BAM file ([bambi decode](https://wtsi-npg.github.io/bambi/#decode))
- Sequencing adapter contamination removal ([biobambam2 bamadapterfind](https://manpages.ubuntu.com/manpages/focal/en/man1/bamadapterfind.1.html))
- BAM to CRAM conversion ([samtools split](http://www.htslib.org/doc/samtools-split.html))
- Data cleaning and QC processes prior to alignment and SNP calling ([cram to bam](https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline/-/blob/develop/workflows/pipeline-subworkflows/cram-to-bam.nf))
- Alignment to full reference / amplicon panel ([realignment](https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline/-/blob/develop/workflows/pipeline-subworkflows/realignment.nf))
- SNP calling using GATK4 tools haplotypecaller and genotypegvcf joint genotyping protocol ([genotyping](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4))

# How to run? (quick and dirty)

Assuming it is running on the farm

**1. Setup**

clone the repository
```
git clone https://gitlab.internal.sanger.ac.uk/malariagen1/ampseq-pipeline.git
cd ./ampseq-pipeine/

```
build the containers
```
cd containers/
bash ./buildContainers.sh
```

WARN: Currently, containers must be built by the user and the pipeline assume the sif files are present at a specific location. Unfortunately, is not possible to build it on the farm, but it can be build on a local machine and copy the .sif files to the farm. The containers will be pull from a registry in the future, but for now that's what we have.  

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
s3_launch_uuid : <str> a universally unique id which will be used to fetch data from s3, if is not provided, the pipeline will not retrieve miseq runs from s3
s3_bucket_input: <str> s3 bucket name to fetch data from

upload_to_s3:<bool> sets if needs to upload output data to an s3 bucket
s3_bucket_output:<str> s3 bucket name to upload data to
```

## Current To Do [1]
- [x] Core replica pipeline
- [x] iRODS
- [x] S3 upload / download
- [x] read counts
- [ ] Genotyping

## Usage

## Support
[who should someone talk to regarding the maintenance and usage of the pipeline]

## Authors and acknowledgment
Show our appreciation to those who have contributed to the project.

## License
[add a license]

## Project status
[prototype]
