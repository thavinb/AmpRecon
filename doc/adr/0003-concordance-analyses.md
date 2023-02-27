# Concordance Analyses

Date: 2023-02-27

## Status
Accepted

## Context
For the ampseq pipeline we decided that an acceptance criteria for completion of a component should be that a concordant analysis should be completed against the inputs/outputs from the Current Production Pipeline (CPP). This means that any given component, provided the same input as the CPP, should produce exactly the same output on a consistent basis.

## Decision
As a team we agreed that an acceptance criteria for any component would be as follows:

### GRC Components
1. CPP genotype files (/lustre/scratch125/gsu/malariagen/production/resources/pipeline_resources/plasmodium/falciparum/amplicon/datasets/genotype_tsvs/) will be used as inputs to be run through any given component. **All genotype files will be run through the component.**
2. Pipeline outputs will be compared to the CPP outputs found in the megaGRC.
    * Care must be taken to ensure that **all** sample IDs found within the megaGRC are found within the pipeline output (e.g. if merging dataframes ensure an outer join is used to ensure records missing from one dataframe are included when datasets are merged, if iterating through datasets ensure that the input set is iterated and that a check is done to ensure a record is presne in the pipeline output).
    * All records should be identical between the two implementations, if they are not they must be justified. Concordance analysis will not be approved unless records with different values are justified/corrected.
        * An example of a justifiable discordance: if the CPP outputs values in a different order to the development pipeline (as occurs in speciation), this could be corrected by sorting the values correctly within your analysis and the sorted values compared.
    * Concordance analyses will maintain the following standards
        * Notebook
            * To the greatest extent possible analyses should be carried out and recorded inside jupyter notebooks
            * Where an operation which is not appropriate for a notebook is carried out, a markdown cell should be present which clearly lays out: the commands that were run (with full paths), the location the command was run in and any other pertinent information
        * Documented
            * A short explainer of the analysis should be present at the top of the notebook along with a short analysis plan
            * The notebook should be divided into sections with markdown header cells according to the analysis plan at the top of the notebook
            * Code should be well commented to ensure clarity in your methodology is conveyed
            * Use of a file from CPP should be documented with a URL to denote where it came from
        * Reproducible
            * Full paths should be provided for all files which are not checked into the repository
            * Prior to submitting for review, restart kernel and re-run notebook from top to bottom (ensures that the same could be done by a reviewer)
        * Robust
            * Analyses should be checked by the developer prior to submitting for review
            * Any potential sources of error (e.g. artifacts upon merging datasets) should be checked (e.g. check no. lines in both dataframes to be merged and the merged dataframe etc)

    **N.B. It has been noted that a potential distinction in this process may be the COI estimation, as this has stochastic elements for which it is unclear whether they could be replicated 100% with the production pipeline. In this event a slightly different approach would likely be taken by comparing distributions of output values. Exact solution is TBD.**

## Consequences
These concordance analyses will provide us with robust eveidence of the effectiveness of our components which can provide us with confidence about their accuracy and provide evidence of their effectiveness which can be presented to stakeholders.