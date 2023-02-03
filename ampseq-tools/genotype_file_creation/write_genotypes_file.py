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

            # Mask - Drop records at particular positions
            non_masked_records  = [record for record in [self._mask_record(record, chromKey_df) for record in vcf_reader] if record != "Masked"]

            # Lift over - Update record co-ordinates
            lifted_over_records = [self._lift_over_record_coordinates(record, chromKey_df) for record in non_masked_records]

            # Filter - Remove low coverage genotypes and adjust depths
            records_genotypes_depths = [self._filter_genotypes(record, vcf_reader) for record in lifted_over_records]

            # Format record
            genotype_file_rows.extend(self._format_record(filtered[0], filtered[1], filtered[2]) for filtered in records_genotypes_depths)

        # Write to genotype file
        genotypes_df = pd.DataFrame(genotype_file_rows, columns = ["Chr", "Loc", "Gen", "Depth", "Filt"])
        genotypes_df.to_csv(self.output_file_name, sep="\t", encoding="utf-8", index=False)

    def _mask_record(self, record, chromKey_df):
        '''
        Labels VCF records at specified positions as "Masked".
        '''
        chromKey_row = self._match_chromKey_row(record, chromKey_df)
        if chromKey_row.iloc[0]["Mask"] == 1:
            return "Masked"
        else:
            return record

    def _lift_over_record_coordinates(self, record, chromKey_df):
        '''
        Updates the chromosome and locus of a supplied VCF record to match those in an associated data frame row.
        '''
        chromKey_row = self._match_chromKey_row(record, chromKey_df)
        record.CHROM = chromKey_row.iloc[0]["Chromosome"]
        record.POS = chromKey_row.iloc[0]["Locus"]
        return record

    def _match_chromKey_row(self, record, chromKey_df):
        '''
        Matches a VCF record to its associated row in the chromKey data frame
        '''
        try:
            matching_chromKey_row = chromKey_df.loc[(chromKey_df["Chrom_ID"] == record.CHROM) & (chromKey_df["VarPos"] == int(record.POS))]
            return matching_chromKey_row
        except IndexError:
            logging.error(f"Failed retrieve matching row from chromKey file for record co-ordinates: {record.CHROM}:{record.POS}")
            exit(1)

    def _filter_genotypes(self, record, vcf_reader):
        '''
        Applies filters to the VCF record genotypes and depths.
        '''
        # If a genotype call exists for this position then continue filtering, otherwise call this position as "missing".
        try:
            # The sample name from this VCF reader
            sample = vcf_reader.samples
            # The genotype call for this record under that sample name in the VCF reader
            call = record.genotype(sample[0])
            # Retrieve genotype call and its associated depths
            genotype = call['GT'] or '.'
            depth = int(call['DP'] or 0)
            allele_depths = call['AD'] or [0, 0, 0, 0]
            allele_depths = allele_depths if isinstance(allele_depths, list) else [allele_depths]
        except:
            return record, "-", "-"
        # Ensure read counts match the alleles

        # Get all alleles for that record - start with the reference
        alleles = [allele for allele in ([record.REF] + record.ALT) if allele != None]
        # Test if we have enough coverage to make a call
        if int(depth) < self.min_total_depth:
            return record, "-", "-"
        if len(alleles) > 1:
            # If there are multiple alleles, filter out alleles with insufficient reads counts or read proportion
            filtered_alleles = []
            filtered_depth = 0
            filtered_allele_depths = []
            for index in range(0, len(alleles)):
                allele_depth = allele_depths[index]
                allele_proportion = allele_depth / depth
                if (allele_depth >= self.het_min_allele_depth) and (allele_proportion >= self.het_min_allele_proportion):
                    filtered_alleles.append(alleles[index])
                    filtered_depth += allele_depth
                    filtered_allele_depths.append(allele_depth)
            # Update allele count and repeat the coverage test
            read_count = 0
            for index in range(0, len(filtered_alleles)):
                read_count += filtered_allele_depths[index]
                print(read_count)
            if int(read_count) < self.min_total_depth:
                return record, "-", "-"
            alleles = filtered_alleles
            allele_depths = filtered_allele_depths
        return record, alleles, allele_depths

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
