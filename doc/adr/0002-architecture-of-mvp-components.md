# 2. Architecture of MVP Components

Date: 2023-01-26

## Status

Accepted

## Context

For the AmpSeq pipeline, we (ARD) need to choose a high-level architecture pattern that helps us replicate the scientific methods of the current production pipeline (https://ssg-confluence.internal.sanger.ac.uk/pages/viewpage.action?pageId=153040844). A modular approach has proven to be an efficient choice in detangling the scientific functionalities of the current production pipeline and logic used to generate GRCs provided to partners. For the MVP, the current scientific rationale, being ported to the MVP, has been documented (https://docs.google.com/document/d/1g_4H4WdoYVYjWanSsYPYBiMq8J1RMWcgcuz8zMdZYOs/edit#heading=h.h6mpjh1o37fg) prior to the start of AmpSeq downstream analysis portion development.
 
Through investigation of several codebases and collaboration with the Parasite analysis team we have identified the following MVP components: <br>
a. Drug Resistance <br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;a1. drug resistance genotyping and haplotype calculation <br> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;a2. kelch13 genotyping <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;a3. plasmepsin <br> 
b. Barcode Generation <br> 
c. Complexity of Infection (RealMcCOIL) <br> 
d. Speciation <br> 
 
We are aiming for (1) a modular pipeline that (2) is not limited to only amplicon sequencing data and will process variant data from both amplicon and WGS runs. To accomplish modularity and validate the scientific method of all components we are developing each component in isolation under the guidelines that all components require the same input file and output is in similar format. Once all components have been reviewed, validated, and agreed upon as a team we will locate shared utility functions and wire components together.

## Decision

As a team, we have agreed to follow the design principles listed below to standardise our development process of the AmpSeq downstream analysis portion and fulfill our goal of a modular pipeline that produces the same data quality as the current production pipeline using the same scientific methods. 
 
1. To fulfill our aim of a pipeline that doesn't care that the data is amplicon sequencing, we will run the pipeline starting from a VCF that comes from a WGS run. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1a. merge per-sample GRC1, GRC2, and Spec VCFs to generate a per sample genotype file that serves as input for all components. <br>
2. For all components, the key for the output data is SampleID. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2a. Output everything in compressed format <br>
3. For input and component outputs, maintain the same naming of the columns. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E.g.: genotype, haplotypes, etc. <br>
4. Minimise file output quantity <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4a. All components are to operate at a level that takes batches as well as per sample. When we begin to do the wiring with Nf we can decide how we want to batch them. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4b. We will have to answer whether serial processing and then reporting as a batch is efficient. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4c. The batching is done on a per MiSeq run basis. Once we are at the phase of wiring process, we can use groovy to automate this batching. <br>
5. The output of each component will be parsed out into 1 final GRC <br>
6. Because the pipeline invovles coding logic, the scientific functionality of each component should be in Python, not groovy. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5a. Think of, and document, what needs to be implemented. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5b. Write a self-contained Python unit of functionality. <br>
7. Upon all components being code reviewed, validated, and agreed upon as a team we will search for shared utility functions and wire components together using groovy (Nextflow). <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7a. Groovy is to be used to glue all components (processes) together. <br>
8. Base image (containers) to be used for development. Everyone starts from the same Python image. <br>
9. For the COI component, the RealMcCOIL will have its own container for all of its dependencies. <br>
10. When writing a Nf process, every process should have an nf-test test. <br>
11. Because our goal is to make this pipeline open-source, Anything Sanger specific is to be modularised. <br>

## Consequences

The design principles we are adhering to help us catch any structuring mistakes early on, minimising risk of technical debt. It allows us to move fast in building the codebase of the first iteration of the MVP. Moving fast in the beginning will allow us to allocate more time to testing and validation. 


