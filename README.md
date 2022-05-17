# AmpSeq-pipeline

The AmpSeq pipeline is a new pipeline developed from scratch to be used by the Amplicon team ([add link]()). [add some background information of the Amplicon project]

## Current Dev Stage: pre_prototype_v0.0

Currently, our main objective is set the piepline structure in a way that is easy to read, understand, maintain, use and alter.
To start the discussions, a mock pipeline was developed and the current implementation just assume the user have text files with named as _<something>.bcl_ and just rename it to a different format.

## Roadmap

1. **mock_pipeline (set the overall structure)**
2. prototype (do real analyses)
3. alpha (to be used by dt5 on isolated actual data analyses scenarios)
4. beta (to be used on a set of real data and stress the code)
5. 1.0 (production ready)

## Current To Do [0.1]

- [x] ~Write basic nextflow framework~
- [ ] add roadmap figure to readme
- [ ] get actual analyses pipeline for the prototype phase
- [x] ~add legacy code Jon codes and modules (mostly for easy access) (DT3)~


## How to deploy?
[add install instructions]

## Usage

**1. Load relevant modules**
```
module load nextflow/21.04.1-5556
```

**2. Running the pipeline**

The pipeline takes as input a directory containing .bcl files. The full path to this directory should be specified in the config file.

Run the pipeline
```
nextflow run main.nf -c pipeline.config
```

## Development Guidelines

Here a set of guidelines for the development process are set just to ensure common understanding between anyone who needs to handle this pipeline in the future.
**THESE ARE NOT RULES** and should be reviewed and changed ~(or even not followed at all)~  if necessary to achieve the main goal of the pipeline.

### Overall Structure
The pipeline is divided into **steps**, defined as any set of interconected **process**, and can be treated as isolated "boxes" which gets some input(s) and generates some output file(s).

 **Steps should be agnostics of what happens to the input file or the output file.**
   Steps should be modular enough to be **changed or removed without any direct effect on the execution of remaining steps** of the workflow.

**Processes** are the actual analyses/command lines that needs to be executed inside "the box".
The processes inside the steps can be aware and/or assume some process needs to be runned before or after it actual execution.

[add real figure of the diagram bellow]

```
        |<- INPUT_PIPELINE
|-- PIPELINE ------------------------------------------------------|
|__ STEPS LAYER __|__ PROCESS LAYER _______________________________|
|- Step1 <- (INPUT_PIPELINE)
|                 |-> proc_1a(step1_input)
|                 |-> proc_1b(output_proc1a)
|                 ...
|                 |-> proc_1N(output_proc1a, output_1i, output_1j)
|-----> Output_Step1
              |
|- Step2 <- (Input)
|                 |-> proc_2a(input_step2)
|                 |-> proc_2b(output_proc2a,input_step2)
|                 ...
|                 |-> proc_2N(output_proc2a, output_2i, ..., output_2j)
|->         Output_Step2
|...
|- StepN <- (Output_StepN-1)
|                 |-> proc_Na(Input)
|                 |-> proc_2b(output_proc2a)
|                 ...
|                 |-> proc_2N(output_proc2a, output_2i, ..., output_2j)
|-------------------------------------------------------------------|
        |-> OUTPUT_PIPELINE
```

### general coding conventions
Set some arbitrary coding conventions is highly advisable to make the code more readable and, therefore, maintainable.
No need to be overly restricted about this, but here are some general conventions for variable and function naming.

Ideally, **variable names** which belongs exclusively to the **step layer** should start with upper case letters (ex: Step1, OutputStep1 or Input_Step4).
vairable names for the **process layer** should only contain lower case letters (ex: run_bwa_index, bamfiles, fastq_channel).

Stuff that run something (functions, process or steps) should have verbs on their names (ex: RunBwaIndex, remap_reads, align_to_ref)

**THIS DEV TEAM BELIEVES IN COMMENTS**. Add comments everywhere, someone in the near future that has close to zero understanding of your code will have to read it.
Most likely, this person will be you.

### Tests
[add guidelines for tests - files and how to run the tests -]


## Support
[who should someone talk to regarding the maintenance and usage of the pipeline]

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
[add a license]

## Project status
[prototype]

## Badges
[add badges]
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.
