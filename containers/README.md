# Containers

This directory holds containers to be used on by the ampseq pipeline.



## How to mount/pull the containers?
Use the bash script (Assumes singularity is available and can be called at the CLI via "_singularity_")

```{bash}
/bin/bash buidContainer.sh
```

PS: this script will **not work at FARM**.
Easiest solution is mount the container outside farm and copy sif files


---

## What are the containers for?

* **bambi**
[Bambi](https://github.com/wtsi-npg/bambi) is set of programs to manipulate SAM/BAM/CRAM files, using HTSLIB.
By default, the container clones the latest version from public github.

Tools currenlty in use: i2b, decode, select.
The processes relying on this container are: **[add list of process]**

* **biobambam2**

[BioBamBam2](https://gitlab.com/german.tischler/biobambam2) is package which contains some tools for processing BAM files.

    bamsormadup: parallel sorting and duplicate marking
    bamcollate2: reads BAM and writes BAM reordered such that alignment or collated by query name
    bammarkduplicates: reads BAM and writes BAM with duplicate alignments marked using the BAM flags field
    bammaskflags: reads BAM and writes BAM while masking (removing) bits from the flags column
    bamrecompress: reads BAM and writes BAM with a defined compression setting. This tool is capable of multi-threading.
    bamsort: reads BAM and writes BAM resorted by coordinates or query name
    bamtofastq: reads BAM and writes FastQ; output can be collated or uncollated by query name

Tools currenlty in use: bamcollate2, bamsort
The processes relying on this container are: **[add list of process]**

* **staden**

[Staden](http://staden.sourceforge.net) is a tool not updated since 2016 and the code is available on sourceforge.
Currently we rely on [scramble](http://staden.sourceforge.net/manual/spin_unix_toc.html#SEC93) at process **[add process name]**
Scramble is a function which produces a version of a given DNA or protein sequence in which the characters are randomly reordered. i.e. the new sequence has the same length and composition as the original but with the characters in a random order.

* **samtools**

[SAMtools](http://www.htslib.org) are two things, 1) a suite of programs for interacting with high-throughput sequencing data. It consists of three separate repositories:
    Samtools: Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
    BCFtools: Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants
    HTSlib: A C library for reading/writing high-throughput sequencing data
This container only contain the Samtools report, not the whole suite.

The processes relying on this container are: **[add list of process]**

* **bwa**

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM.

The processes relying on this container are: **[add list of process]**
