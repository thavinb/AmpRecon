# COI component

Here a script to handle input and output for the REALMcCOIL and a container recipe for the whole component is provided

The `grc_process_mccoil_io.py` can:

1) write input files for the REAL McCOIL from sample barcodes `tsv` and a barcode definition `json` file. Both should follow the same format as set for the Barcode component of the `grc_tools`;
2) write `coi.grc` from the `_summary.txt` output of the REAL McCOIL.

PS: Inputs and outputs of the original The REAL McCOIL should work just fine, however, here we use a [custom version of The REAL McCOIL](https://github.com/AMarinhoSN/THEREALMcCOIL). Small tweeks were necessary to containarize this software and improve usability by adding a command line interface was necessary.

## USAGE

```
usage: grc_process_mccoil_io.py [-h] [-write_mccoil_in] [-write_coi_grc] [--barcodes_files BARCODES_FILES [BARCODES_FILES ...]] [--barcode_def_file BARCODE_DEF_FILE]
                                [--mccoil_sum_file MCCOIL_SUM_FILE [MCCOIL_SUM_FILE ...]] [--output_file OUTPUT_FILE]

A script to write McCOIL input from barcodes tsv files and write grc files from McCOIL output

optional arguments:
  -h, --help            show this help message and exit
  -write_mccoil_in      prepare McCOIL input files from barcodes
  -write_coi_grc        prepare McCOIL input files from barcodes
  --barcodes_files BARCODES_FILES [BARCODES_FILES ...]
                        Path(s) to barcode tsv file for a batch of samples
  --barcode_def_file BARCODE_DEF_FILE
                        Path to a json file with barcodes definitions
  --mccoil_sum_file MCCOIL_SUM_FILE [MCCOIL_SUM_FILE ...]
                        Path to McCOIL summary output file
  --output_file OUTPUT_FILE
                        Path for the McCOIl input or grc file to be written (default: <./McCOIL_in.tsv | ./coi.grc >). If more than one input file is provided, the output will have the same
                        name of the respective input file added of a suffix ("_mccoil_in.tsv" | "_coi.grc")
```

A common usage of this component, assuming it is running The REALMcCOIL inside the container, would look like:

```
# write McCOIL in
python3 grc_process_mccoil_io.py -write_mccoil_in \
            --barcodes_files $BARCODES_IN --config $BARCODE_DEF \
            --output_file $RUN_ID.tsv

# run McCOIL
Rscript /app/THEREALMcCOIL/runMcCOIL.R -i ${RUN_ID}.tsv \
            --totalRun ${NTOTAL} --totalBurnIn ${NBURN} --seed 123456 \
            --outPrefix $RUN_ID --maxCOI 20 --M0 5

# write coi grc
python3 grc_process_mccoil_io.py -write_coi_grc \
            --mccoil_sum_file ${RUN_ID}_summary.txt \
            --output_file ${RUN_ID}_out.grc
 - `WriteCOIgrc.py`: write a grc containing sample and coi results based on output summary file from McCOIL
```
