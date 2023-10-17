# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- **[Feature]**: Added fastq entry point.
- **[Feature]**: NF-test based unit test of the workflow `reads_to_variants.nf`.
- **[Feature]**: Added new "run_locally" profile for pipeline execution without containers.

### Changed
- **[Change]**: Updated README.
- **[Removed]**: Removal of post processing read counts step.
- **[Improvement]**: Checked validity of a custom path provided via `mccoil_repopath`.


## [1.1.1] - 2023-10-02

### Added
- **[Feature]**: use ampseq resources files provided by submodule by default.
- **[Feature]**: config settings specific for Pf and Pv.
- **[Feature]**: NF-test based unit test of the workflow `miseq_to_reads.nf`.
- **[Feature]**: NF-test based unit test of the workflow `sanger_irods_to_reads.nf`.
- **[Feature]**: NF-test based unit test of the workflow `variants_to_grcs.nf`.
- **[Feature]**: Switches allow kelch-13 and plasmepsin GRC creation steps to be switched off.

### Changed
- **[improvement]**: ENA cram files are now published on the output dir.
- **[improvement]**: Output dir now is structured on subdirs.
- **[improvement]**: Added missing parameter check calls.
- **[improvement]**: Allow kelch-reference file to not be specified if `no_kelch` flag set.
- **[improvement]**: Changed the pipeline to only output a single GRC: this GRC contains all the output data produced by variants_to_grcs.

### Fixed
- **[Bug]**: bug for runids with "_T" affecting read counts file fixed.
- **[Bug]**: upload data to s3 function is fixed
- **[Bug]**: Column names for Pf were hard-coded into metadata addition script, this blocked the usage of the pipeline with Pv.
  - **Fix**: Column ordering was abstracted to GRC settings file and metadata addition script was modified to include config file and read ordering from there.

### Removed
- **[Change]**: Removed s3 input processes.
- **[Change]**: Removed requirement for RunID to be an integer
- **[Change]**: Removed GRC formatting and column sorting


## [1.0.1] - 2023-07-13

### Fixed

- **[Bug]** the manifest for in-country cannot have empty/NA values on any columns. This is a problem for test runs, which will always have empty values for metadata columns.

  - **Fix**: The `validate_manifest.py` was changed and now checks for empty/NA values only on the required columns for the pipeline to run properly (`'sample_id'`, `'primer_panel'`, `'barcode_number'`,`'barcode_sequence'`). Any extra column will not be checked.

- **[Bug]** post processing of read couts files breaks when using test run ids standarad naming (`{run_id}_T`).
  - **Fix**: The `post_proc_read_counts.py` now removes any `_T` before recovering the lanelets.

## [1.0.0] - 2023-06-15

- This version of the pipeline is functionally equivalent to the Current Production Pipeline, but represent massive improvements on various front, which include (and is not limited to) code base maintainability, flexibility, scalability and portability.
