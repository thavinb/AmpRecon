
## Stage 1 pipeline deployment

### Introduction
This document details what is required to run the Stage-1 ASPIS pipeline. More information about what the Stage-1 ASPIS pipeline is can be found [here](readme.md)

### Required software
The following packages are required to run the Stage-1 pipeline

- [Bambi](http://wtsi-npg.github.io/bambi/)
- [biobambam2](https://github.com/gt1/biobambam2/tree/master/src/programs)
- [bio-bwa](http://bio-bwa.sourceforge.net/)
- [SamTools](http://www.htslib.org/doc/samtools.html)
- [BCF Tools](https://samtools.github.io/bcftools/bcftools.html)
- [Staden - sCRAMble](https://github.com/jkbonfield/io_lib)

### Hardware requirements

The main hardware requirement for the Stage-1 pipeline is hard-drive space. The output of the sequencing machine are ~62 GB. Despite each step in the pipeline reducing the size of data it stays in the GBs until the halfway point so 100 GB can be quickly used up for a single job / batch.

I have not ran a full scale test but so far I am able to work through the steps of the pipeline at a reasonable rate using a m1.medium virtual machine hosted on OpenStack. This machine has 4 VCPUs and ~33 GB of ram

As all the OpenStack flavours have hard drive size that increase inline with memory a volume was used to provide the additional storage space required

### OpenStack installation guide

For this deployment I heavily rely on the packages and recipes created by wtsi-npg, hosted on [this repo](https://github.com/wtsi-npg/npg_conda)

1. Create a new instance on OpenStack using a m1.medium virtual machine running Ubuntu 18
2. Install miniconda
 * `mkdir miniconda`
 * `cd miniconda/``
 * `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`
 * `bash Miniconda3-latest-Linux-x86_64.sh`
3. Check out npg_conda
 * `cd ~`
 * `git clone https://github.com/wtsi-npg/npg_conda.git`
4. Follow npg_conda [HOWTO](https://github.com/wtsi-npg/npg_conda/blob/devel/doc/HOWTO.md)
 * `nano ~/.condarc` Enter the following
```
 show_channel_urls: true
channels:
  - https://dnap.cog.sanger.ac.uk/npg/conda/prod/Ubuntu/18.04/
  - defaults
```
 * 'conda install conda-build'
5. Create conda environment and install required packages
 * `conda create --name pipe-tools`
 * `conda activate pipe-tools`
 * `conda install bambi`
 * `conda install biobambam2`
 * `conda install bwa-aligner`
 * `conda install samtools`
 * `conda install bcftools`
6. Mount additional storage - Mounting instructions can be found [here](https://access.redhat.com/documentation/en-US/Red_Hat_Enterprise_Linux_OpenStack_Platform/2/html/Getting_Started_Guide/ch16s03.html)
 * Create a new volume on OpenStack (or reuse an existing one)
 * Attach it to the Stage-1 Pipeline instance
 * SSH into the instance
 * Find the disk id `ls /dev/disk/by-id`
 * If it is a new volume create a filesystem `mkfs.ext4 [disk path]`
 * Create a folder to mount the drive to `mkdir -p /mnt/bcl`
 * Mount the drive `mount [disk path] /mnt/bcl`
