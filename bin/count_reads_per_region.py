#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

import argparse
import sys
import os
import re
import csv
import logging
from collections import defaultdict
from typing import List, Tuple
import pysam

# Set up logging to file
logging.basicConfig(
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


class Counters:
    """
    A dataclass for storing count values
    """

    def __init__(
        self,
        zero: int = 0,
        ten: int = 0,
        ten_f: int = 0,
        ten_r: int = 0,
        ten_frag: int = 0,
        ten_fragboth: int = 0,
    ):
        self.zero = zero
        self.ten = ten
        self.ten_f = ten_f
        self.ten_r = ten_r
        self.ten_frag = ten_frag
        self.ten_fragboth = ten_fragboth


class QC:
    """Calculate QC counts using pysam. Output the results as a csv to stdout"""

    DEFAULT_QUALITY_THRESHOLD = 10
    _DESIGN_FILE_REGEX = re.compile(r"^([\w.]+):(\d+)-(\d+)$")
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
        self._out_handle = None
        self._output_hdr = (
            "Rpt,Region,Amplicon_name,Total_reads,Total_region_reads,Region_reads,Perc of total reads,"
            "Perc of mapped to region reads,Total_region_reads MQ>={qual},Region_reads MQ>={qual},"
            "Region_reads 1 MQ>={qual},Region_reads 2 MQ>={qual},Perc Region_reads 1 MQ>={qual},"
            "Perc of mapped to region reads MQ>={qual},Region_fragments represented MQ>={default_qual},"
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
        if self._out_handle and self._out_handle != sys.stdout:
            self._out_handle.close()

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

    def _count_reads_in_region(self, bam_file: pysam.AlignmentFile, region: str, 
                               mapq_threshold: int = 0, region_start: int = None, 
                               region_end: int = None) -> Tuple[int, int, int, int, int]:
        """
        Count reads in a specific region using pysam
        
        Returns:
            - total_reads: Total number of reads in region
            - forward_reads: Forward strand reads (flag 0x40)  
            - reverse_reads: Reverse strand reads (flag 0x80)
            - fragments: Number of unique fragment names
            - fragments_both: Number of fragments with both pairs
        """
        total_reads = 0
        forward_reads = 0
        reverse_reads = 0
        fragment_counts = defaultdict(int)
        
        try:
            # Parse region if it's in chr:start-end format
            if ':' in region:
                chrom, coords = region.split(':', 1)
                if '-' in coords:
                    start, end = map(int, coords.split('-'))
                else:
                    start, end = None, None
            else:
                chrom = region
                start, end = None, None
            
            logger.debug(f"Fetching reads from {chrom}:{start}-{end} with MAPQ>={mapq_threshold}")
            
            # Iterate through reads in the region
            # Use separate parameters for pysam.fetch() 
            # Note: pysam uses 0-based coordinates, so convert from 1-based design file
            pysam_start = (start - 1) if start else None
            pysam_end = end if end else None
            
            for read in bam_file.fetch(contig=chrom, start=pysam_start, stop=pysam_end):
                # Skip reads based on flags (equivalent to -F0xB04)
                if (read.is_unmapped or read.is_secondary or read.is_supplementary):
                    continue
                
                # Apply MAPQ filter
                if read.mapping_quality < mapq_threshold:
                    continue
                
                # Apply position filtering if region coordinates are provided and region_start is set
                if region_start is not None and region_end is not None:
                    read_start = read.reference_start + 1  # pysam is 0-based, SAM is 1-based
                    
                    # Use template length for calculating read end (matching original logic)
                    # NOTE: Do NOT use abs() - template_length can be negative for reverse read
                    template_length = read.template_length if read.template_length else 0
                    if template_length != 0:
                        read_end = read_start + template_length - 1
                    else:
                        read_end = read.reference_end if read.reference_end else read_start + len(read.query_sequence)
                    
                    # Check if read is within 5bp of region boundaries
                    x = read_start - region_start
                    y = read_end - region_end
                    
                    # Must have mate on same reference (RNEXT == "=") and be within boundaries
                    # Original script checks fields[6] == "=" which is "mate on same reference"
                    # In pysam: read.next_reference_id == read.reference_id means mate is on same reference
                    mate_on_same_ref = (read.next_reference_id == read.reference_id)
                    if not (mate_on_same_ref and (x <= 5 and x >= -5) and (y <= 5 and y >= -5)):
                        continue
                    
                    total_reads += 2  # Count both pairs for position-filtered regions
                else:
                    # For regions without specific coordinates, count each read
                    total_reads += 1
                    
                    # Count forward/reverse reads (only for non-position-filtered)
                    if read.is_read1:
                        forward_reads += 1
                    elif read.is_read2:
                        reverse_reads += 1
                    
                # Count fragments for all cases
                fragment_counts[read.query_name] += 1
                        
        except Exception as e:
            logger.warning(f"Error processing region {region}: {e}")
            return 0, 0, 0, 0, 0
        
        # Calculate fragment statistics
        fragments = len(fragment_counts)
        fragments_both = sum(1 for count in fragment_counts.values() if count == 2)
        
        logger.debug(f"Region {region}: total_reads={total_reads}, fragments={fragments}, fragments_both={fragments_both}")
        
        return total_reads, forward_reads, reverse_reads, fragments, fragments_both

    def _get_total_reads(self, bam_file: pysam.AlignmentFile) -> int:
        """Get total number of reads (excluding only secondary, supplementary - includes unmapped!)"""
        total = 0
        # Use until_eof=True to get all reads including unmapped, then filter
        logger.debug("Counting reads manually with until_eof=True to match samtools -F0xB00...")
        for read in bam_file.fetch(until_eof=True):
            # Exclude only secondary and supplementary (equivalent to -F0xB00)
            if not (read.is_secondary or read.is_supplementary):
                total += 1
        logger.debug(f"Total reads from manual count with until_eof (equiv to -F0xB00): {total}")
        return total

    def run(self):
        """Run the process"""
        logger.info("Starting QC with pysam")

        # Read input files
        region_config = self._read_text_file(self.design_file)
        plex_list = self._read_csv_file(self.plex_file)

        # Write qc file header
        self._out_handle.write(self._output_hdr)

        for plex_entry in plex_list:
            logger.debug("Processing sample: {}".format(plex_entry[0]))

            rpt = plex_entry[0]
            bam_file_path = os.path.join(self.input_dir, rpt + ".bam")
            
            if not os.path.exists(bam_file_path):
                logger.warning(f"BAM file {bam_file_path} not found, skipping...")
                continue

            # Open BAM file
            try:
                bam_file = pysam.AlignmentFile(bam_file_path, "rb")
                logger.debug(f"Successfully opened BAM file: {bam_file_path}")
                
                # Check if BAM file has index
                try:
                    bam_file.check_index()
                    logger.debug("BAM file has index")
                except:
                    logger.warning(f"BAM file {bam_file_path} may not be indexed")
                    
            except Exception as e:
                logger.error(f"Failed to open BAM file {bam_file_path}: {e}")
                continue

            total_reads = self._get_total_reads(bam_file)
            total_region_reads = 0
            total_region_reads_mq = 0
            total_fragments = 0
            count = defaultdict(Counters)

            logger.info(f"Processing {len(region_config)} regions for sample {rpt}")

            for line in region_config:
                parts = line.split(",")
                region = parts[0]
                amplicon_name = parts[1] if len(parts) > 1 else region
                
                # Check if region has coordinates
                m = self._DESIGN_FILE_REGEX.match(region)
                if m:
                    region_start = int(m.group(2))
                    region_end = int(m.group(3))
                    use_position_filter = True
                else:
                    region_start = None
                    region_end = None
                    use_position_filter = False

                logger.debug(f"Processing region: {region} (amplicon: {amplicon_name})")

                # Count reads without MAPQ filter
                num_all, fwd_all, rev_all, frag_all, frag_both_all = self._count_reads_in_region(
                    bam_file, region, 0, region_start, region_end
                )

                # Count reads with MAPQ filter
                num_mq, fwd_mq, rev_mq, frag_mq, frag_both_mq = self._count_reads_in_region(
                    bam_file, region, self.mapq, region_start, region_end
                )
                
                logger.debug(f"Region {region}: all_reads={num_all}, mq_reads={num_mq}, position_filter={use_position_filter}")

                # For position-filtered regions, adjust the fragment counts
                if use_position_filter:
                    fwd_mq = num_mq // 2
                    rev_mq = num_mq // 2
                    frag_mq = num_mq // 2
                    frag_both_mq = num_mq // 2

                total_region_reads += num_all
                total_region_reads_mq += num_mq
                total_fragments += frag_both_mq

                count[region] = Counters(
                    zero=num_all,
                    ten=num_mq,
                    ten_f=fwd_mq,
                    ten_r=rev_mq,
                    ten_frag=frag_mq,
                    ten_fragboth=frag_both_mq,
                )

            # Calculate unmapped/unused reads
            unmap_reads = total_reads - total_region_reads
            unmap_reads_mq = total_reads - total_region_reads_mq

            count[self._UNMAP_OR_UNUSED_KEY] = Counters(
                zero=unmap_reads, ten=unmap_reads_mq
            )

            logger.info(
                f"Total for {bam_file_path} => {total_region_reads} ({total_reads}) - MQ>={self.mapq} {total_region_reads_mq}"
            )

            # Write results
            for line in region_config:
                parts = line.split(",")
                r = parts[0]
                amplicon_name = parts[1] if len(parts) > 1 else r
                
                c = count[r]
                self._out_handle.write(
                    ",".join(
                        str(field)
                        for field in [
                            rpt,  # Rpt
                            r,  # Region
                            amplicon_name,  # Amplicon_name
                            total_reads,  # Total_reads
                            total_region_reads,  # Total_region_reads
                            c.zero,  # Region_reads
                            "0" if total_reads < 1 else f"{(c.zero / total_reads) * 100:.2f}",  # Perc of total reads
                            "0" if total_region_reads < 1 else f"{(c.zero / total_region_reads) * 100:.2f}",  # Perc of mapped to region reads
                            total_region_reads_mq,  # Total_region_reads MQ>=threshold
                            c.ten,  # Region_reads MQ>=threshold
                            c.ten_f,  # Region_reads 1 MQ>=threshold
                            c.ten_r,  # Region_reads 2 MQ>=threshold
                            "0" if c.ten < 1 else f"{(c.ten_f / c.ten) * 100:.2f}",  # Perc Region_reads 1 MQ>=threshold
                            "0" if total_region_reads_mq < 1 else f"{(c.ten / total_region_reads_mq) * 100:.2f}",  # Perc of mapped to region reads MQ>=threshold
                            c.ten_frag if c.ten_frag else "0",  # Region_fragments represented MQ>=threshold
                            c.ten_fragboth if c.ten_fragboth else "0",  # Region_fragments both MQ>=threshold
                            f"{(c.ten_fragboth / ((total_reads / 2) + 0.00000001)) * 100:.2f}",  # Perc of total fragments
                            f"{(c.ten_fragboth / total_fragments) * 100:.2f}" if total_fragments > 0 else "0",  # Perc of mapped to region fragments
                        ]
                    )
                    + "\n"
                )

            # Write unmapped/unused line
            uc = count[self._UNMAP_OR_UNUSED_KEY]
            self._out_handle.write(
                ",".join(
                    str(field)
                    for field in [
                        rpt,  # Rpt
                        self._UNMAP_OR_UNUSED_KEY,  # Region
                        self._UNMAP_OR_UNUSED_KEY,  # Amplicon_name
                        total_reads,  # Total_reads
                        total_region_reads,  # Total_region_reads
                        uc.zero,  # Region_reads
                        "0" if total_reads < 1 else f"{(uc.zero / total_reads) * 100:.2f}",  # Perc of total reads
                        "",  # Perc of mapped to region reads
                        total_region_reads_mq,  # Total_region_reads MQ>=threshold
                        uc.ten,  # Region_reads MQ>=threshold
                        "",  # Region_reads 1 MQ>=threshold
                        "",  # Region_reads 2 MQ>=threshold
                        "",  # Perc Region_reads 1 MQ>=threshold
                        "",  # Perc of mapped to region reads MQ>=threshold
                        "",  # Region_fragments represented MQ>=threshold
                        "",  # Region_fragments both MQ>=threshold
                        f"{(((total_reads / 2) - total_fragments) / ((total_reads / 2) + 0.00000001)) * 100:.2f}",  # Perc of total fragments
                        "",  # Perc of mapped to region fragments
                    ]
                )
                + "\n"
            )

            bam_file.close()

        logger.info("Finished QC")


def get_arguments():
    """Read program arguments"""
    parser = argparse.ArgumentParser(
        description="Run QC on a set of BAM files using pysam. Logging is written to a qc.log file."
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
        help="Directory containing the sorted/indexed BAM files.",
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
            f"ERROR: design file {args.design_file} does not exist or is zero length\n"
        )
        return error_status
    if not os.path.exists(args.input_dir):
        sys.stderr.write(
            f"ERROR: input_dir directory {args.input_dir} does not exist\n"
        )
        return error_status

    # Log to file
    fh = logging.FileHandler(args.log_file, mode="w")
    formatter = logging.Formatter(
        "[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Run QC
    try:
        with QC(
            args.design_file, args.plex_file, args.input_dir, args.mapq, args.output_file
        ) as qc_runner:
            qc_runner.run()
    except Exception as e:
        logger.error(f"QC failed: {e}")
        return error_status

    return 0


if __name__ == "__main__":
    sys.exit(main())