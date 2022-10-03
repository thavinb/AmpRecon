# AmpSeq-pipeline

The AmpSeq pipeline is a new pipeline developed from scratch to be used by the Amplicon team ([add link]()). [add some background information of the Amplicon project]

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
download_from_s3: <bool> sets if needs to download data from a S3 bucket (in-country entry point)
uuid: <str> a universally unique id which will be used to fetch data from s3
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
