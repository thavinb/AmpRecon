#!/usr/bin/env python3
import argparse
import sys
import os
import re
import csv
import time
import subprocess
import collections
import logging
from multiprocessing.pool import ThreadPool
from os import path
from typing import List

# Set up logging to file
logging.basicConfig(
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


class Counters:
    """
    A dataclass for storing count values

    Can't use the new dataclass annotation due to the Python version being used
    """

    def __init__(
        self,
        zero: float = 0.0,
        ten: float = 0.0,
        ten_f: float = 0.0,
        ten_r: float = 0.0,
        ten_frag: float = 0.0,
        ten_fragboth: float = 0.0,
    ):
        self.zero = zero
        self.ten = ten
        self.ten_f = ten_f
        self.ten_r = ten_r
        self.ten_frag = ten_frag
        self.ten_fragboth = ten_fragboth


class QC:
    """Calculate QC counts on a given set of cram files. Output the results as a csv to stdout"""

    DEFAULT_QUALITY_THRESHOLD = 10

    _CRAM_FILE_NAME_FORMAT = "{}"
    _DESIGN_FILE_REGEX = re.compile(r"^([\w.]+):(\d+)-(\d+)$")
    _CMD_BASE = "samtools view "
    _UNMAP_OR_UNUSED_KEY = "UNMAP_or_UNUSED"

    def __init__(
        self,
        design_file: str,
        plex_file: str,
        input_dir: str,
        mapq: int,
        output_file: str,
    ) -> None:
        self.design_file = design_file
        self.plex_file = plex_file
        self.input_dir = input_dir
        self.mapq = mapq
        self.output_file = output_file
        self._pool = ThreadPool(processes=5)
        self._out_handle = None
        self._output_hdr = (
            "Rpt,Region,Total_reads,Total_region_reads,Region_reads,Perc of total reads,"
            "Perc of mapped to region reads,,Total_region_reads MQ>={qual},Region_reads MQ>={qual},"
            "Region_reads 1 MQ>={qual},Region_reads 2 MQ>={qual},Perc Region_reads 1 MQ>={qual},"
            "Perc of mapped to region reads MQ>={qual},,Region_fragments represented MQ>={default_qual},"
            "Region_fragments both MQ>={qual},Perc of total fragments,Perc of mapped to region fragments\n".format(
                qual=mapq, default_qual=self.DEFAULT_QUALITY_THRESHOLD
            )
        )

    def __enter__(self):
        """Makes the class function like a context manager"""
        self._out_handle = (
            self.output_file and open(self.output_file, "w") or sys.stdout
        )
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the context manager"""
        self.cleanup()

    def cleanup(self):
        """Close the output file and thread pool"""
        if self._out_handle:
            self._out_handle.close()
        if self._pool:
            self._pool.close()

    @staticmethod
    def _read_text_file(file_name: str) -> List[str]:
        """Read in a text file"""
        lines = []
        with open(file_name) as fh:
            for line in fh:
                li = line.strip()
                # Ignore any comment lines
                if not li.startswith("#"):
                    lines.append(li)
        return lines

    @staticmethod
    def _read_csv_file(file_name: str) -> List[List[str]]:
        """Read in a given comma separated file"""
        lines = []
        with open(file_name) as f:
            csv_file = csv.reader(f, delimiter=",")
            for line in csv_file:
                lines.append(line)
        return lines

    @classmethod
    def _run_samtools(cls, args: str) -> List[str]:
        """Run samtools and wait for an exit status. Return the output"""
        attemps_left = 5
        while attemps_left>0:
            try:
                process = subprocess.run(
                    cls._CMD_BASE + args,
                    shell=True,
                    universal_newlines=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )

                stdout_lines = process.stdout.splitlines()
                stderr_lines = process.stderr
                if len(stderr_lines):
                    logger.warning(stderr_lines)

                return stdout_lines

            except subprocess.CalledProcessError as err:
                attemps_left -= 1 ;
                if attemps_left>0:
                    time.sleep(5.0) # Sleep 5 seconds
                    continue
                logger.exception(err)
                raise RuntimeError(
                    "[ {} ] returned error: {} and a status of {}".format(
                        cls._CMD_BASE + args,
                        str(err).replace("\n", " "),
                        str(err.returncode),
                    )
                ) from err

    def _start_task(self, args: List[str]):
        """Execute an async thread in the pool"""
        arg_string = " ".join(args)
        logger.debug("Executing task: " + self._CMD_BASE + arg_string)
        return self._pool.apply_async(
            self._run_samtools, (arg_string,), error_callback=self._handle_error
        )

    @staticmethod
    def _calculate_samtools_count(
        text: List[str], region_start: str, region_end: str
    ) -> int:
        """Take a list of text output lines from samtools and calculate the count"""
        count = 0
        region_start_int = int(region_start)
        region_end_int = int(region_end)
        for line in text:
            fields = line.split("\t")
            x = int(fields[3]) - region_start_int
            y = ((int(fields[3]) + int(fields[8])) - 1) - region_end_int
            if (x <= 5 and x >= -5) and (fields[6] == "=") and (y <= 5 and y >= -5):
                count += 2

        return count

    @staticmethod
    def _get_samtools_count(text: List[str]) -> int:
        """Take a list of text output lines from samtools and get the integer count from it"""
        return 0 if not text else int(text[0])

    @staticmethod
    def _calculate_samtools_frag_counts(text: List[str]) -> (int, int):
        """Take a list of text output lines from samtools and calculate fragment counts from it"""
        fragments = {}
        y = 0
        for line in text:
            fields = line.split("\t")
            frag_id = fields[0]
            if frag_id in fragments:
                fragments[frag_id] += 1
            else:
                fragments[frag_id] = 1
        x = len(fragments)
        for key, value in fragments.items():
            if fragments[key] == 2:
                y += 1
        return x, y

    def _handle_error(self, e: Exception):
        """Handle fatal errors from samtools"""
        logger.fatal(str(e))
        self.cleanup()
        # Make sure we exit threads...
        os._exit(-1)

    def run(self):
        """Run the process"""
        logger.info("Starting QC")

        # Read input files
        region_config = self._read_text_file(self.design_file)
        plex_list = self._read_csv_file(self.plex_file)

        # Write qc file header
        self._out_handle.write(self._output_hdr)

        for plex_entry in plex_list:

            logger.debug("Next region : {}".format(plex_entry))

            rpt = self._CRAM_FILE_NAME_FORMAT.format(
                plex_entry[0]
            )
            cram_file = path.join(self.input_dir, rpt + ".bam")

            total0 = 0
            total10 = 0
            frc = 0.0
            count = collections.defaultdict()

            for region in region_config:
                m = self._DESIGN_FILE_REGEX.match(region)

                if m:
                    region_start = m.group(2)
                    region_end = m.group(3)
                else:
                    region_start = None
                    region_end = None

                # Calculate qc counts
                if region_start is not None:
                    task1 = self._start_task(["-F0xB04", cram_file, region])
                    task2 = self._start_task(
                        ["-q", str(self.mapq), "-F0xB04", cram_file, region]
                    )
                    # Get results from threads
                    num = self._calculate_samtools_count(
                        task1.get(), region_start, region_end
                    )
                    num2 = self._calculate_samtools_count(
                        task2.get(), region_start, region_end
                    )
                    num3 = num2 / 2
                    num4 = num2 / 2
                    frag = num2 / 2
                    frag_both = num2 / 2

                else:
                    task1 = self._start_task(["-c", "-F0xB04", cram_file, region])
                    task2 = self._start_task(
                        ["-c", "-q", str(self.mapq), "-F0xB04", cram_file, region]
                    )
                    task3 = self._start_task(
                        [
                            "-c",
                            "-q",
                            str(self.mapq),
                            "-F0xB04",
                            "-f0x0040",
                            cram_file,
                            region,
                        ]
                    )
                    task4 = self._start_task(
                        [
                            "-c",
                            "-q",
                            str(self.mapq),
                            "-F0xB04",
                            "-f0x0080",
                            cram_file,
                            region,
                        ]
                    )
                    task5 = self._start_task(
                        ["-q", str(self.mapq), "-F0xB04", cram_file, region]
                    )
                    # Get results from threads
                    num = self._get_samtools_count(task1.get())
                    num2 = self._get_samtools_count(task2.get())
                    num3 = self._get_samtools_count(task3.get())
                    num4 = self._get_samtools_count(task4.get())
                    frag, frag_both = self._calculate_samtools_frag_counts(task5.get())

                total0 += num
                total10 += num2
                count[region] = Counters(
                    zero=num,
                    ten=num2,
                    ten_f=num3,
                    ten_r=num4,
                    ten_frag=frag,
                    ten_fragboth=frag_both,
                )
                frc += frag_both

            task = self._start_task(["-c", "-F0xB00", cram_file])
            numa = self._get_samtools_count(task.get())
            count[self._UNMAP_OR_UNUSED_KEY] = Counters(
                zero=(numa - total0), ten=(numa - total10)
            )

            logger.info(
                "Total for {} => {} ({}) - MQ>={} {}".format(
                    cram_file, total0, numa, self.mapq, total10
                )
            )

            # Write results to stdout

            for r in region_config:
                self._out_handle.write(
                    ",".join(
                        str(field)
                        for field in [
                            rpt,
                            r,
                            numa,
                            total0,
                            int(count[r].zero),
                            "0"
                            if numa < 1
                            else "{:.2f}".format((count[r].zero / numa) * 100),
                            "0"
                            if total0 < 1
                            else "{:.2f}".format((count[r].zero / total0) * 100),
                            "",
                            total10,
                            int(count[r].ten),
                            int(count[r].ten_f),
                            int(count[r].ten_r),
                            "0"
                            if count[r].ten < 1
                            else "{:.2f}".format((count[r].ten_f / count[r].ten) * 100),
                            "0"
                            if total10 < 1
                            else "{:.2f}".format((count[r].ten / total10) * 100),
                            "",
                            int(count[r].ten_frag) if count[r].ten_frag else "0",
                            int(count[r].ten_fragboth)
                            if count[r].ten_fragboth
                            else "0",
                            "{:.2f}".format(
                                (count[r].ten_fragboth / ((numa / 2) + 0.00000001))
                                * 100
                            ),
                            "{:.2f}".format((count[r].ten_fragboth / frc) * 100)
                            if frc > 0
                            else "0",
                        ]
                    )
                    + "\n"
                )

            self._out_handle.write(
                ",".join(
                    str(field)
                    for field in [
                        rpt,
                        self._UNMAP_OR_UNUSED_KEY,
                        numa,
                        total0,
                        int(count[self._UNMAP_OR_UNUSED_KEY].zero),
                        "0"
                        if numa < 1
                        else "{:.2f}".format(
                            (count[self._UNMAP_OR_UNUSED_KEY].zero / numa) * 100
                        ),
                        "",
                        "",
                        total10,
                        int(count[self._UNMAP_OR_UNUSED_KEY].ten),
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "{:.2f}".format(
                            (((numa / 2) - frc) / ((numa / 2) + 0.00000001)) * 100
                        ),
                        "",
                    ]
                )
                + "\n"
            )

        logger.info("Finished QC")


def get_arguments():
    """Read program arguments"""
    parser = argparse.ArgumentParser(
        description="Run QC on a set of CRAM files. Logging is written to a qc.log file."
    )
    parser.add_argument(
        "--design_file",
        "-d",
        dest="design_file",
        type=str,
        required=True,
        help="The input QC configuration file.",
    )
    parser.add_argument(
        "--plex_file",
        "-p",
        dest="plex_file",
        type=str,
        required=True,
        help="Input QC plex file.",
    )
    parser.add_argument(
        "--input_dir",
        "-i",
        dest="input_dir",
        type=str,
        required=True,
        help="Directory containing the sorted/indexed CRAM files.",
    )
    parser.add_argument(
        "--mapq",
        "-q",
        dest="mapq",
        type=int,
        required=False,
        default=QC.DEFAULT_QUALITY_THRESHOLD,
        help="Use a non-default quality threshold. Default: 10",
    )
    parser.add_argument(
        "--output",
        "-o",
        dest="output_file",
        type=str,
        required=False,
        help="An optional output file name. Default: use stdout",
    )
    parser.add_argument(
        "--log_file",
        "-l",
        dest="log_file",
        type=str,
        required=False,
        default="qc.log",
        help="An optional log file name. Default: qc.log",
    )
    return parser


def main():
    # Parse program arguments
    args = get_arguments().parse_args()
    error_status = 255

    # Check arguments
    if not os.path.exists(args.design_file) or os.path.getsize(args.design_file) == 0:
        sys.stderr.write(
            "ERROR: design file {} does not exist or is zero length\n".format(
                args.design_file
            )
        )
        return error_status
    if not os.path.exists(args.plex_file) or os.path.getsize(args.plex_file) == 0:
        sys.stderr.write(
            "ERROR: plex file {} does not exist or is zero length\n".format(
                args.plex_file
            )
        )
        return error_status
    if not os.path.exists(args.input_dir):
        sys.stderr.write(
            "ERROR: input_dir directory {} does not exist\n".format(str(args.input_dir))
        )
        return error_status

    # Log to file
    fh = logging.FileHandler(args.log_file, mode="w")
    formatter = logging.Formatter("[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Run QC
    with QC(
        args.design_file, args.plex_file, args.input_dir, args.mapq, args.output_file
    ) as qc_runner:
        qc_runner.run()


if __name__ == "__main__":
    sys.exit(main())


