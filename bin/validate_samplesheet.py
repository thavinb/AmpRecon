#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

"""
Manifest validation script.

This script validates a manifest file for Nextflow pipelines. It performs several checks,
including:
- Verifying the presence of required headers.
- Ensuring file paths have the correct format (e.g., .cram or .fastq).
- Checking for duplicate entries based on an internal pipeline ID.
- Validating primer panel names against a provided list.

The script supports two execution modes: 'irods' and 'fastq', each with its own
specific validation rules.
"""

import logging
from typing import List, Tuple
import pandas as pd
from sys import argv, exit

# Configure the logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

def load_dataframes(
    mnf_flpath: str,
    panel_set_flpath: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load the manifest and panel dataframes.

    Args:
        mnf_flpath: Path to the manifest file.
        panel_set_flpath: Path to the panel set file.

    Returns:
        A tuple containing the manifest DataFrame and the panel set DataFrame.
    """

    logging.info("Loading Manifest and Panel Setting Files")
    try:
        df = pd.read_csv(mnf_flpath, sep='\t')
        pnl_df = pd.read_csv(panel_set_flpath, sep=',')
        return df, pnl_df
    except FileNotFoundError as e:
        logging.error(f"A required file was not found - {e.filename}")
        exit(1)
    except pd.errors.ParserError as e:
        logging.error(f"Failed to parse a file - {e}")
        exit(1)

def get_config(execution_mode: str) -> Tuple[List[str], Tuple, List[str]]:
    """
    Set configuration based on execution mode.

    Args:
        execution_mode: The mode of execution ('cram' or 'fastq').

    Returns:
        A tuple containing:
        - A list of required headers.
        - A tuple of required file formats.
        - A list of columns containing file paths.
    """

    if execution_mode == 'cram':
        logging.info(f"Execution Mode: {execution_mode}")
        return ["sample_id", "primer_panel", "cram_path"], (".cram",), ["cram_path"]
    elif execution_mode == 'fastq':
        logging.info(f"Execution Mode: {execution_mode}")
        return ["sample_id", "primer_panel", "fastq1_path", "fastq2_path"], ('.fastq', '.fq', '.fastq.gz', '.fq.gz'), ["fastq1_path", "fastq2_path"]
    else:
        logging.error(f"Unsupported execution mode '{execution_mode}'. Choose 'cram' or 'fastq'.")
        exit(1)

def validate_headers(
    df: pd.DataFrame,
    required_headers: List[str]
) -> bool:
    """
    Validate that all required headers are present in the DataFrame.

    Args:
        df: The manifest DataFrame.
        required_headers: A list of expected header names.

    Returns:
        True if all headers are present, False otherwise.
    """
    logging.info("Validate Column Headers")

    columns_lst = df.columns.tolist()
    missing_headers = [h for h in required_headers if h not in columns_lst]

    if missing_headers:
        logging.error(f"The following expected headers were not found: {missing_headers}")
        return False
    return True

def validate_file_formats(
    df: pd.DataFrame,
    name_path_columns: List[str],
    required_format: Tuple
) -> bool:
    """
    Validate that file paths in specified columns have the required format.

    Args:
        df: The manifest DataFrame.
        name_path_columns: A list of columns to check for file paths.
        required_format: A tuple of acceptable file extensions.

    Returns:
        True if all file paths have the correct format, False otherwise.
    """
    logging.info(f"Validate File Extensions: {required_format}" )

    invalid_ext_list = set([])
    for col in name_path_columns:
        invalid_mask = ~df[col].apply(lambda x: x.endswith(required_format))
        if invalid_mask.any():
            invalid_idx = df[invalid_mask].index
            invalid_ext_list = invalid_ext_list.union(set(invalid_idx))

    if invalid_ext_list:
        logging.error(f"{len(invalid_ext_list)} lines were found with invalid extentions in the manifest (i.e., {required_format}): {invalid_ext_list}")
        return False

    return True

def validate_uniqueness(
    df: pd.DataFrame,
    execution_mode: str
) -> bool:
    """
    Validate that internal pipeline IDs are unique.

    Args:
        df: The manifest DataFrame.
        execution_mode: The mode of execution ('irods' or 'fastq').

    Returns:
        True if all IDs are unique, False otherwise.
    """
    logging.info("Ensure non-duplicated records")
    logging.info("fastq: sampleid_panelid; cram: filename_sampleid_panelid")

    def gen_pipeline_internal_id(row: pd.Series) -> str:
        """Create a unique ID for pipeline usage."""
        if execution_mode == 'cram':
            # Handle potential non-string or NaN values
            cram_path = str(row["cram_path"]) if pd.notna(row["cram_path"]) else ""
            simple_name = cram_path.split("/")[-1].split(".")[0]
            return f"{simple_name}_{row['sample_id']}_{row['primer_panel']}"
        elif execution_mode == 'fastq':
            return f"{row['sample_id']}_{row['primer_panel']}"
        return ""

    df["internal_pipeline_id"] = df.apply(gen_pipeline_internal_id, axis=1)

    duplicated_bool = df["internal_pipeline_id"].duplicated(keep=False)
    duplicated_pipe_ids = df[duplicated_bool]

    if not duplicated_pipe_ids.empty:
        logging.error(f"{len(duplicated_pipe_ids)} lines in the manifest were found duplicated: {duplicated_pipe_ids.index}")
        return False

    return True

def validate_primer_panels(
    df: pd.DataFrame,
    required_panels: List[str]
) -> bool:
    """
    Validate that primer panel values are valid.

    Args:
        df: The manifest DataFrame.
        required_panels: A list of valid primer panel names.

    Returns:
        True if all panel names are valid, False otherwise.
    """
    logging.info(f"Validate manifest's primer panel name against the panel setting file: {required_panels}")

    invalid_ppnls = df[~df["primer_panel"].isin(required_panels)]

    if not invalid_ppnls.empty:
        logging.error(f"{len(invalid_ppnls)} lines in the manifest have primer panel names that are not in the provided panel settings file: {invalid_ppnls.index.to_list()}, ")
        logging.error(f"Invalid provided primer_panel: {invalid_ppnls.primer_panel.unique()}")
        return False

    return True

def main():
    if len(argv) != 4:
        print("Usage: validate_manifest.py <manifest_file> <panel_set_file> <execution_mode>")
        exit(1)

    mnf_flpath, panel_set_flpath, execution_mode = argv[1], argv[2], argv[3]
    errors_found = 0

    df, pnl_df = load_dataframes(mnf_flpath, panel_set_flpath)

    required_headers, required_format, name_path_columns = get_config(execution_mode)
    required_panels = pnl_df["panel_name"].to_list()

    # Run all validation checks
    if not validate_headers(df.copy(), required_headers):
        errors_found += 1

    if not validate_file_formats(df.copy(), name_path_columns, required_format):
        errors_found += 1

    if not validate_uniqueness(df.copy(), execution_mode):
        errors_found += 1

    if not validate_primer_panels(df.copy(), required_panels):
        errors_found += 1

    if errors_found > 0:
        logging.error(f"{errors_found} errors found, please check the manifest and/or panel setting files")
        exit(1)
    logging.info(":: DONE ::")

if __name__ == "__main__":
    main()
