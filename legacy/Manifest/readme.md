# Stage 1 Manifest - V0.1.2

## Introduction
The stage 1 manifest contains all the information required to run both stage 1 and 2 of the ASPIS pipelines. It consists of two parts, a header section and a data section. It is formatted as a comma separated value file (CSV) for easy access and creation.

## Header
This section holds information that is relevant to the pipeline and to every sample within the job. Information in this section are stored as key value pairs across two columns per key value pair. The start of this section is denoted by `[Header]` in the first column on a line on its own

### Values
- **Version** This should be 0.1.2. Any other value will be rejected by the Pipeline.
- **Workflow** This should be the workflow entered into the [MiSeq](https://www.illumina.com/systems/sequencing-platforms/miseq.html) via the SampleSheet
- **Pipeline** This holds the id of the stage 2 pipeline that you with to process the data with. For example `nihr-pf`
- **Pipeline Version** This holds the version of the stage 2 pipeline. For example for `nihr-pf` the current version is `0.1.0`
- **Date** Date of submission in YYYY-MM-dd format. For example `2019-05-03`
- **Lane** Lane used during sequencing

## Data
The section holds information per sample to be processed. It is denoted by `[Data]` in the first column on a line on its own. On the next line the following column headers should be added `LIMS ID, SIMS ID, Index, Ref, Barcode Sequence`. Each sample should be entered on a separate line under this header.

### Values
- **LIMS ID** You id for a sample. Must be unique
- **SIMS ID** A Roma ID. Must be unique
- **Index** Samples index
- **Ref** - Reference used for sample (PFA_GRC1_v1.0, PFA_Spec, PFA_GRC2_v1.0)
- **Barcode Sequence** - Tag A and Tag B concatenated with a `-` hyphen (ATCACGTT-TATAGCCT)
- **Well** - Holds the well reference
- **Plate** - Holds the plate id

## Example
```
[Header]
Version,	0.1.2
Workflow,	LibraryQC
Pipeline, nihr_pf
Pipeline Version, 0.1.0
Date,	2019-06-06
Lane,	1
[Data]
LIMS ID,	SIMS ID,	Index,	Ref, 	Barcode Sequence, Well, Plate
3429STDY7977738,RCN15176,1,PFA_GRC1_v1.0,ATCACGTT-TATAGCCT,A1,DN572422E
3429STDY7977739,RCN15177,2,PFA_GRC1_v1.0,CGATGTTT-TATAGCCT,B1,DN572422E
3429STDY7977740,RCN15178,3,PFA_GRC1_v1.0,TTAGGCAT-TATAGCCT,C1,DN572422E
3429STDY7977741,RCN15179,4,PFA_GRC1_v1.0,TGACCACT-TATAGCCT,D1,DN572422E
3429STDY7977742,RCN15180,5,PFA_GRC1_v1.0,ACAGTGGT-TATAGCCT,E1,DN572422E
3429STDY7977743,RCN15181,6,PFA_GRC1_v1.0,GCCAATGT-TATAGCCT,F1,DN572422E
3429STDY7977744,RCN15182,7,PFA_GRC1_v1.0,CAGATCTG-TATAGCCT,G1,DN572422E
3429STDY7977745,RCN15183,8,PFA_GRC1_v1.0,ACTTGATG-TATAGCCT,H1,DN572422E
3429STDY7977746,RCN15184,9,PFA_GRC1_v1.0,GATCAGCG-TATAGCCT,A2,DN572422E
...
```
## Tools
In the tool section of manifest there are a set of manifest focused tools

- `sanger2manifest.py` - Creates a manifest using data from Sequence Scape. Requires a Run id
- `manifest2taglist.py` - Accepts a manifest and creates a taglist required by Stage 1
- `manifest2manifest.py` - Accespts a manifest and uses ROMA to create a Stage 2 manifest

## Validation
A validation tool does not currently exist but you can use `manifestParser.py` to confirm that a manifest is valid.
