#!/usr/bin/env python3

import argparse
import pandas as pd
import logging
import gzip
import vcf

class GenotypeFileWriter:
    '''
    Merges a list of VCF files into a single genotype .tsv file. Updates the co-ordinates of the rows. Drops rows at specified loci.
    '''
    def __init__(self, vcf_file_path_list, out_file_name, chromKey_file_path, min_total_depth, het_min_allele_depth, het_min_allele_proportion):
        self.vcf_list = vcf_file_path_list
        self.output_file_name = out_file_name
        self.chromKey_file = chromKey_file_path
        self.min_total_depth = float(min_total_depth)
        self.het_min_allele_depth = float(het_min_allele_depth)
        self.het_min_allele_proportion = float(het_min_allele_proportion)

    def write_genotype_file(self):
        '''
        Merges supplied VCF files into a single genotype .tsv file.
        Updates VCF record co-ordinates to match those found within ChromKey data frame.
        SNPs to be masked are not written to the resulting file.
        '''
        chromKey_df = pd.read_csv(self.chromKey_file, sep="\t", header=0)
        genotype_file_rows = []
        for vcf_file in self.vcf_list:
            vcf_reader = vcf.Reader(filename=vcf_file)

            # Match - Get row associated with record from chromKey file
            records_and_chromKey_rows = [self._match_chromKey_row(record, chromKey_df) for record in vcf_reader]

            # Mask - Drop records at particular positions
            non_masked_data  = [record for record in [self._mask_record(row[0], row[1]) for row in records_and_chromKey_rows] if record != "Masked"]

            # Lift over - Update record co-ordinates
            lifted_over_records = [self._lift_over_record_coordinates(row[0], row[1]) for row in non_masked_data]

            # Filter - Remove low coverage genotypes and adjust depths
            records_genotypes_depths = [self._filter_genotypes(record, vcf_reader) for record in lifted_over_records]

            # Format record
            genotype_file_rows.extend(self._format_record(filtered[0], filtered[1], filtered[2]) for filtered in records_genotypes_depths)

        genotypes_df = pd.DataFrame(genotype_file_rows, columns = ["Chr", "Loc", "Gen", "Depth", "Filt"])
        genotypes_df.to_csv(self.output_file_name, sep="\t", encoding="utf-8", index=False)

    def _match_chromKey_row(self, record, chromKey_df):
        '''
        Matches a VCF record to its associated row in the chromKey data frame
        '''
        try:
            matching_chromKey_row = chromKey_df.loc[(chromKey_df["Chrom_ID"] == record.CHROM) & (chromKey_df["VarPos"] == int(record.POS))]
            return record, matching_chromKey_row
        except IndexError:
            logging.error(f"Failed retrieve matching row from chromKey file for record co-ordinates: {record.CHROM}:{record.POS}")
            exit(1)

    def _mask_record(self, record, chromKey_row):
        '''

        Labels VCF records at specified positions as "Masked".
        '''
        if chromKey_row.iloc[0]["Mask"] == 1:
            return "Masked"
        else:
            return record, chromKey_row

    def _lift_over_record_coordinates(self, record, chromKey_row):
        '''
        Updates the chromosome and locus of a supplied VCF record to match those in an associated data frame row.
        '''
        record.CHROM = chromKey_row.iloc[0]["Chromosome"]
        record.POS = chromKey_row.iloc[0]["Locus"]
        return record

    def _filter_genotypes(self, record, vcf_reader):
        '''
        Returns only the VCF record alleles with enough coverage to make a call.
        '''
        # If a genotype call and its depths exists for this position then continue filtering, otherwise call position as "missing"
        try:
            sample = vcf_reader.samples
            call = record.genotype(sample[0])

            # Retrieve genotype and depths for this record
            genotype = call['GT'] or '.'
            depth = int(call['DP'] or 0)
            allele_depths = call['AD'] or [0, 0, 0, 0]
            allele_depths = allele_depths if isinstance(allele_depths, list) else [allele_depths]

            # Get all alleles for this record - start with the reference
            alleles = [allele for allele in ([record.REF] + record.ALT) if allele != None]
        except:
            return record, "-", "-"

        # Sanity Check: Ensure the number of allele depths match the number of alleles
        if len(allele_depths) != len(alleles):
            logging.error(f"Error parsing VCF {sample[0]}: at position {record.CHROM}:{record.POS} AD field contains {len(allele_depths)} values for {len(alleles)} alleles.")
            exit(1)

        # Test that this record has enough coverage to make a call
        if int(depth) < self.min_total_depth:
            return record, "-", "-"

        # If there are multiple alleles, filter out alleles with insufficient depths and or read proportions
        if len(alleles) > 1:
            alleles, allele_depths, updated_depth = self._filter_individual_alleles(alleles, allele_depths, depth)
            # Test whether remaining alleles have enough coverage to make a call
            if int(updated_depth) < self.min_total_depth:
                return record, "-", "-"

        return record, alleles, allele_depths

    def _filter_individual_alleles(self, alleles, allele_depths, depth):
        '''
        Return alleles and depths that exceed depth and proportion thresholds
        '''
        filtered_alleles = []
        filtered_allele_depths = []
        updated_depth = 0
        # Iterate over each the alleles of this record
        for index in range(0, len(alleles)):
            # Retrieve allele depth and calculate what proportion of reads at this position it comprises
            allele_depth = allele_depths[index]
            allele_proportion = allele_depth / depth
            if (allele_depth >= self.het_min_allele_depth) and (allele_proportion >= self.het_min_allele_proportion):
                # Retain alleles with depths that exceed threshold and which comprise more than minimum proportion of all the reads at this position
                filtered_alleles.append(alleles[index])
                filtered_allele_depths.append(allele_depth)
                updated_depth += allele_depth
        return filtered_alleles, filtered_allele_depths, updated_depth

    def _format_record(self, record, genotypes, allele_depths):
        '''
        Creates a row for output to a genotype file.
        '''
        genotypes_formatted = ",".join(str(genotype) for genotype in genotypes)
        allele_depths_formatted = ",".join(str(depth) for depth in allele_depths) if isinstance(allele_depths, list) else allele_depths
        filter_value = record.FILTER[0] if len(record.FILTER) != 0 else "PASS"
        return [record.CHROM, str(record.POS), str(genotypes_formatted), str(allele_depths_formatted), str(filter_value)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf_list", "-i", help="List of VCF files to merge, mask SNPs in and update co-ordinates for.", required=True, type=str, nargs="+", default=[])
    parser.add_argument("--output_file_name", "-o", help="Base name (no extension) for the output genotypes file.", required=True)
    parser.add_argument("--chromKey_file", "-m", help="File containing SNPs to mask in merged VCF.", required=True)
    parser.add_argument("--min_total_depth", "-c", help=".", required=True)
    parser.add_argument("--het_min_allele_depth", "-a", help=".", required=True)
    parser.add_argument("--het_min_allele_proportion", "-p", help=".", required=True)
    args = parser.parse_args()
    liftover_genotypes_write_file = GenotypeFileWriter(args.input_vcf_list, args.output_file_name, args.chromKey_file, args.min_total_depth, args.het_min_allele_depth, args.het_min_allele_proportion)
    liftover_genotypes_write_file.write_genotype_file()
