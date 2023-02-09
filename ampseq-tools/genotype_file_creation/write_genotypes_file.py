#!/usr/bin/env python3

import argparse
import pandas as pd
import logging
import vcf

class GenotypeFileWriter:
    '''
    Merges a list of VCF files into a single genotype .tsv file.
    Drops records at specified loci, updates the co-ordinates of the rows and filters out alleles with low coverage.
    '''
    def __init__(self, vcf_file_path_list, out_file_name, chromKey_file_path, chromosome_column, locus_column, min_total_depth, het_min_allele_depth, het_min_allele_proportion):
        self.vcf_list = vcf_file_path_list
        self.output_file_name = out_file_name
        self.chromKey_file = chromKey_file_path
        self.chromosome_column = chromosome_column
        self.locus_column = locus_column
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

        # Iterate over the supplied VCF files
        for vcf_file in self.vcf_list:
            vcf_reader = vcf.Reader(filename=vcf_file)

            # Match - Get row associated with record from chromKey file
            records_matched_to_chromKey_rows = [self._match_chromKey_row(record, chromKey_df) for record in vcf_reader]

            # Mask - Drop records at particular positions
            remaining_records_and_rows  = [record for record in [self._mask_record(row[0], row[1]) for row in records_matched_to_chromKey_rows] if record != "Masked"]

            # Lift over - Update record co-ordinates
            lifted_over_records = [self._lift_over_record_coordinates(row[0], row[1]) for row in remaining_records_and_rows]

            # Filter - Remove low coverage genotypes and adjust depths
            genotypes_depths = [self._filter_genotypes(record[0], vcf_reader) for record in lifted_over_records]

            # Format record
            genotype_file_rows.extend(self._format_record(record[1], record[0], value[1], value[2]) for record, value in zip(lifted_over_records, genotypes_depths))

        # Write formatted records to genotype file
        genotypes_df = pd.DataFrame(genotype_file_rows, columns = ["Amplicon", "Pos", "Chr", "Loc", "Gen", "Depth", "Filt"])
        genotypes_df.to_csv(self.output_file_name, sep="\t", encoding="utf-8", index=False)

    def _match_chromKey_row(self, record, chromKey_df):
        '''
        Matches a VCF record to its associated row in the chromKey data frame
        '''
        try:
            matching_chromKey_row = chromKey_df.loc[(chromKey_df[str(self.chromosome_column)] == record.CHROM) & (chromKey_df[str(self.locus_column)] == int(record.POS))]
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
        amplicon = [chromKey_row.iloc[0]["Chrom_ID"], chromKey_row.iloc[0]["VarPos"]] 
        record.CHROM = chromKey_row.iloc[0]["Chromosome"]
        record.POS = chromKey_row.iloc[0]["Locus"]
        return record, amplicon

    def _filter_genotypes(self, record, vcf_reader):
        '''
        Returns only the VCF record alleles with enough coverage to make a call.
        '''
        # If a genotype call and its depths exists for this position then continue filtering, otherwise call position as "missing"
        try:
            sample = vcf_reader.samples
            call = record.genotype(sample[0])

            # Retrieve genotype and depths for this record
            genotype = call['GT']
            depth = int(call['DP'])
            allele_depths = call['AD']
            allele_depths = allele_depths if isinstance(allele_depths, list) else [allele_depths]

            # Get all alleles for this record - start with the reference
            alleles = [allele for allele in ([record.REF] + record.ALT) if allele != None]
        except:
            return "-", "-"

        # Sanity Check: Ensure the number of allele depths match the number of alleles
        if len(allele_depths) != len(alleles):
            logging.error(f"{call} Error parsing VCF {sample[0]}: at position {record.CHROM}:{record.POS} AD field contains {len(allele_depths)} values for {len(alleles)} alleles.")
            exit(1)

        # Test that this record has enough coverage to make a call
        if int(depth) < self.min_total_depth:
            return "-", "-"

        # If there are multiple alleles, filter out alleles with insufficient depths and read proportions
        if len(alleles) > 1:
            alleles, allele_depths, updated_depth = self._filter_individual_alleles(alleles, allele_depths, depth)
            # Test whether remaining alleles have enough coverage to make a call
            if int(updated_depth) < self.min_total_depth:
                return "-", "-"

        return alleles, allele_depths

    def _filter_individual_alleles(self, alleles, allele_depths, depth):
        '''
        Return alleles and depths that exceed depth and proportion thresholds
        '''
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
            if (allele_depth >= self.het_min_allele_depth) and (allele_proportion >= self.het_min_allele_proportion):
                filtered_alleles.append(alleles[index])
                filtered_allele_depths.append(allele_depth)
                updated_depth += allele_depth

        return filtered_alleles, filtered_allele_depths, updated_depth

    def _format_record(self, amplicon, record, genotypes, allele_depths):
        '''
        Creates a row for output to a genotype file.
        '''
        genotypes_formatted = ",".join(str(genotype) for genotype in genotypes)
        allele_depths_formatted = ",".join(str(depth) for depth in allele_depths) if isinstance(allele_depths, list) else allele_depths
        filter_value = ";".join(record.FILTER) if record.FILTER is not None and len(record.FILTER) > 0 else "PASS"
        return [amplicon[0], amplicon[1], record.CHROM, str(record.POS), str(genotypes_formatted), str(allele_depths_formatted), str(filter_value)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf_list", "-i", help="List of VCF files to merge, mask SNPs in and update co-ordinates for.", required=True, type=str, nargs="+", default=[])
    parser.add_argument("--output_file_name", "-o", help="Base name (no extension) for the output genotypes file.", required=True)
    parser.add_argument("--chromKey_file", "-m", help="File containing SNPs to mask in merged VCF.", required=True)
    parser.add_argument("--chromosome_column_name", "-c", help="Name of the column in the chromKey file to try to match chromosome to.", required=True)
    parser.add_argument("--locus_column_name", "-l", help="Name of the column in the chromKey file to try to match position to.", required=True)
    parser.add_argument("--min_total_depth", "-d", help="The number of reads at a record must exceed this value.", required=True)
    parser.add_argument("--het_min_allele_depth", "-a", help="The number of reads for a particular allele at a heterozygous record must exceed this value.", required=True)
    parser.add_argument("--het_min_allele_proportion", "-p", help="The proportion of reads for particular allele at a heterozygous record must exceed this value.", required=True)
    args = parser.parse_args()
    liftover_genotypes_write_file = GenotypeFileWriter(args.input_vcf_list, args.output_file_name, args.chromKey_file, args.chromosome_column_name, args.locus_column_name, args.min_total_depth, args.het_min_allele_depth, args.het_min_allele_proportion)
    liftover_genotypes_write_file.write_genotype_file()
