# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- **[Feature]**: None

### Changed
- **[improvement]**: ENA cram files are now published on the output dir.
- **[improvement]**: Output dir now is structured on subdirs.
- **[improvement]**: Added missing parameter check calls.
### Fixed

- **[Bug]** the manifest for in-country cannot have empty/NA values on any columns. This is a problem for test runs, which will always have empty values for metadata columns.

  - **Fix**: The `validate_manifest.py` was changed and now checks for empty/NA values only on the required columns for the pipeline to run properly (`'sample_id'`, `'primer_panel'`, `'barcode_number'`,`'barcode_sequence'`). Any extra column will not be checked.

- **[Bug]** post processing of read couts files breaks when using test run ids standarad naming (`{run_id}_T`).
  - **Fix**: The `post_proc_read_counts.py` now removes any `_T` before recovering the lanelets.

## [1.0.0] - 2023-06-15

- This version of the pipeline is functionally equivalent to the Current Production Pipeline, but represent massive improvements on various front, which include (and is not limited to) code base maintainability, flexibility, scalability and portability.