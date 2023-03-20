#!/usr/bin/env python3

import argparse
import logging
import vcf
import csv


class GenotypeFileWriter:
    """
    Merges a list of VCF files into a single genotype .tsv file.
    Drops records at specified loci, updates the co-ordinates of the rows and filters out alleles with low coverage.
    """

    def __init__(
        self,
        vcf_list,
        out_file_name,
        sample_id,
        chromKey_file_path,
        chromosome_column,
        locus_column,
        min_total_depth,
        het_min_allele_depth,
        het_min_allele_proportion,
    ):
        self.vcf_list = vcf_list
        self.output_file_name = out_file_name
        self.sample_id = sample_id
        self.chromKey_file = chromKey_file_path
        self.chromosome_column = chromosome_column
        self.locus_column = locus_column
        self.min_total_depth = float(min_total_depth)
        self.het_min_allele_depth = float(het_min_allele_depth)
        self.het_min_allele_proportion = float(het_min_allele_proportion)

    def write_genotype_file(self):
        """
        Merges supplied VCF files into a single genotype .tsv file.
        Updates VCF record co-ordinates to match those found within ChromKey dictionary.
        SNPs to be masked are not written to the resulting file.
        """
        # Read ChromKey file rows as dictionaries
        chromKey_dict = {}
        chromKey_file = open(self.chromKey_file, newline="")
        for row in csv.DictReader(chromKey_file, delimiter="\t"):
            key = f"{row[str(self.chromosome_column)]}:{row[str(self.locus_column)]}"
            chromKey_dict[key] = dict(row)
        chromKey_file.close()

        formatted_vcf_rows = []

        # Iterate over the supplied VCF files and their records
        for vcf_file in self.vcf_list:

            vcf_reader = vcf.Reader(filename=vcf_file)
            sample = vcf_reader.samples

            for record in vcf_reader:

                # Match - Get row associated with record from chromKey file
                chromKey_row = self._match_chromKey_row(record, chromKey_dict)

                # Mask - Ignore this record if its supposed to be masked
                if chromKey_row.get("Mask") == "1":
                    continue

                # Lift over - Update record co-ordinates
                record, amplicon = self._lift_over_record_coordinates(
                    record, chromKey_row
                )

                # Filter - Remove low coverage genotypes and adjust depths
                genotype, depth = self._filter_genotypes(record, sample)

                # Format record
                formatted_vcf_rows.append(
                    self._format_record(amplicon, record, genotype, depth)
                )

        # Store formatted rows from the VCFs under a key of the chromosome:locus
        genotype_rows_dict = {f"{row[0]}:{row[1]}": row for row in formatted_vcf_rows}

        # Remove positions from chromKey that are to be masked
        filtered_chromKey_dict = dict(
            filter(lambda value: value[1].get("Mask") != "1", chromKey_dict.items())
        )

        # Write output genotype file
        with open(self.output_file_name, "w") as output_genotype_file:
            file_writer = csv.DictWriter(
                output_genotype_file,
                delimiter="\t",
                fieldnames=[
                    "ID",
                    "Amplicon",
                    "Pos",
                    "Chr",
                    "Loc",
                    "Gen",
                    "Depth",
                    "Filt",
                ],
            )
            file_writer.writeheader()

            # Format each chromKey row and add any VCF record data at that position
            for key, row in filtered_chromKey_dict.items():
                formatted_row = self._format_chromKey_row(key, row, genotype_rows_dict)

                # Write formatted chromKey / VCF rows to genotype file
                file_writer.writerow(formatted_row)

    def _match_chromKey_row(self, record, chromKey_dict):
        """
        Matches a VCF record to its associated dictionary in the chromKey dictionary.
        Matches by supplied columns (default is the amplicon and amplicon position)
        """
        try:
            matching_chromKey_row = chromKey_dict.get(f"{record.CHROM}:{record.POS}")
            return matching_chromKey_row
        except IndexError:
            logging.error(
                f"Failed retrieve matching row from chromKey file for record co-ordinates: {record.CHROM}:{record.POS}"
            )
            exit(1)

    def _lift_over_record_coordinates(self, record, chromKey_row):
        """
        Updates the chromosome and locus of a supplied VCF record to match those in an associated dictionary.
        Also returns amplicon and amplicon position.
        """
        amplicon = [chromKey_row.get("Chrom_ID"), chromKey_row.get("VarPos")]
        record.CHROM = chromKey_row.get("Chromosome")
        record.POS = chromKey_row.get("Locus")
        return record, amplicon

    def _filter_genotypes(self, record, sample):
        """
        Returns only the VCF record alleles with enough coverage to make a call.
        """
        # If a genotype call and its depths exists for this position then continue filtering, otherwise call position as "missing"
        try:
            call = record.genotype(sample[0])
            # Retrieve genotype and depths for this record
            depth = int(call["DP"])
            allele_depths = call["AD"]
            allele_depths = (
                allele_depths if isinstance(allele_depths, list) else [allele_depths]
            )

            # Get all alleles for this record - start with the reference
            alleles = [
                allele for allele in ([record.REF] + record.ALT) if allele != None
            ]
        except:
            return "-", "-"

        # Sanity Check: Ensure the number of allele depths match the number of alleles
        if len(allele_depths) != len(alleles):
            logging.error(
                f"{call} Error parsing VCF {sample[0]}: at position {record.CHROM}:{record.POS} AD field contains {len(allele_depths)} values for {len(alleles)} alleles."
            )
            exit(1)

        # Test that this record has enough coverage to make a call
        if int(depth) < self.min_total_depth:
            return "-", "-"

        # If there are multiple alleles, filter out alleles with insufficient depths and read proportions
        if len(alleles) > 1:
            alleles, allele_depths, updated_depth = self._filter_individual_alleles(
                alleles, allele_depths, depth
            )
            # Test whether remaining alleles have enough coverage to make a call
            if int(updated_depth) < self.min_total_depth:
                return "-", "-"

        return alleles, allele_depths

    def _filter_individual_alleles(self, alleles, allele_depths, depth):
        """
        Return alleles and depths that exceed depth and proportion thresholds.
        """
        filtered_alleles = []
        filtered_allele_depths = []
        updated_depth = 0

        # Iterate over the alleles of this record
        for index in range(0, len(alleles)):

            # Retrieve allele depth
            # Calculate what proportion of the total depth at this position this allele's depth comprises
            allele_depth = allele_depths[index]
            allele_proportion = allele_depth / depth

            # Retain alleles with depths that both:
            # 1) Exceed minimum depth threshold
            # 2) Comprise more than minimum proportion of all the reads at this position
            if (allele_depth >= self.het_min_allele_depth) and (
                allele_proportion >= self.het_min_allele_proportion
            ):
                filtered_alleles.append(alleles[index])
                filtered_allele_depths.append(allele_depth)
                updated_depth += allele_depth

        return filtered_alleles, filtered_allele_depths, updated_depth

    def _format_record(self, amplicon, record, genotypes, allele_depths):
        """
        Creates a row for output to a genotype file.
        """
        genotypes_formatted = ",".join(str(genotype) for genotype in genotypes)
        allele_depths_formatted = (
            ",".join(str(depth) for depth in allele_depths)
            if isinstance(allele_depths, list)
            else allele_depths
        )
        filter_value = (
            ";".join(record.FILTER)
            if isinstance(record.FILTER, list) and len(record.FILTER) > 0
            else ("." if record.FILTER == None else "PASS")
        )
        return [
            amplicon[0],
            amplicon[1],
            record.CHROM,
            str(record.POS),
            str(genotypes_formatted),
            str(allele_depths_formatted),
            str(filter_value),
        ]

    def _format_chromKey_row(self, key, row, genotypes_dict):
        """
        Format a chromKey row and adds any VCF record data at that position
        """

        # If a genotype row taken from a VCF matches then use its values, otherwise use "-"
        genotypes_row = genotypes_dict.get(key)
        if genotypes_row != None:
            genotype = genotypes_row[4]
            depth = genotypes_row[5]
            filter = genotypes_row[6]
        else:
            genotype = "-"
            depth = "-"
            filter = "-"

        # Retrieve location data from this chromKey row
        amplicon = row.get("Chrom_ID")
        amplicon_pos = row.get("VarPos")
        chromosome = row.get("Chromosome")
        chr_loc = row.get("Locus")

        row = {
            "ID": str(self.sample_id),
            "Amplicon": amplicon,
            "Pos": str(amplicon_pos),
            "Chr": chromosome,
            "Loc": str(chr_loc),
            "Gen": genotype,
            "Depth": depth,
            "Filt": filter,
        }
        return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vcf_files",
        "-i",
        help="List of VCF files to merge, mask SNPs in and update co-ordinates for.",
        required=True,
        type=str,
        nargs="+",
        default=[],
    )
    parser.add_argument(
        "--output_file",
        "-o",
        help="Name of the output genotypes file.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--sample_id",
        "-s",
        help="Sample ID of the VCF files - output to genotype file.",
        required=True,
    )
    parser.add_argument(
        "--chromKey_file",
        "-m",
        help="File containing SNPs to mask in merged VCF.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--chromosome_column_name",
        "-c",
        help="Name of the column in the chromKey file to try to match chromosome to.",
        type=str,
        default="Chrom_ID",
    )
    parser.add_argument(
        "--locus_column_name",
        "-l",
        help="Name of the column in the chromKey file to try to match position to.",
        type=str,
        default="VarPos",
    )
    parser.add_argument(
        "--min_total_depth",
        "-d",
        help="The number of reads at a record must exceed this value.",
        type=float,
        default=10,
    )
    parser.add_argument(
        "--het_min_allele_depth",
        "-a",
        help="The number of reads for a particular allele at a heterozygous record must exceed this value.",
        type=float,
        default=5,
    )
    parser.add_argument(
        "--het_min_allele_proportion",
        "-p",
        help="The proportion of reads for particular allele at a heterozygous record must exceed this value.",
        type=float,
        default=0.1,
    )
    args = parser.parse_args()
    liftover_genotypes_write_file = GenotypeFileWriter(
        args.vcf_files,
        args.output_file,
        args.sample_id,
        args.chromKey_file,
        args.chromosome_column_name,
        args.locus_column_name,
        args.min_total_depth,
        args.het_min_allele_depth,
        args.het_min_allele_proportion,
    )
    liftover_genotypes_write_file.write_genotype_file()
