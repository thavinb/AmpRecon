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
    def __init__(self, vcf_file_path_list, out_file_name, chromKey_file_path, MinCallReads, MinAlleleReads, MinAlleleProp):
        self.vcf_list = vcf_file_path_list
        self.output_file_name = out_file_name
        self.chromKey_file = chromKey_file_path
        self.min_call_reads = float(MinCallReads)
        self.min_allele_reads = float(MinAlleleReads)
        self.min_allele_proportion = float(MinAlleleProp)

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

            # Filter
            records_genotypes_depths = [self._filter_record(record, vcf_reader) for record in lifted_over_records]

            # Format record
            genotype_file_rows.extend(self._format_record(filtered[0], filtered[1], filtered[2]) for filtered in records_genotypes_depths)

        # Write to genotype file
        genotypes_df = pd.DataFrame(genotype_file_rows, columns = ["Chr", "Loc", "Gen", "Depth", "Filt"])
        genotypes_df.to_csv(self.output_file_name, sep="\t", encoding="utf-8", index=False)
#I may alter this beginning part slightly
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
            chromKey_row = chromKey_df.loc[(chromKey_df["Chrom_ID"] == record.CHROM) & (chromKey_df["VarPos"] == int(record.POS))]
            return chromKey_row
        except IndexError:
            logging.error(f"Failed retrieve matching row from chromKey file for record co-ordinates: {record.CHROM}:{record.POS}")
            exit(1)

    def _filter_record(self, record, vcf_reader):
        '''
        Applies filters to the VCF record genotypes and depths.
        '''
        # If a genotype call exists for this position then continue filtering, otherwise call this position as "missing".
        try:
            # The sample name from this VCF
            sample = vcf_reader.samples
            # The genotype call for this record under that sample name in the VCF reader
            call = record.genotype(sample[0])
            # Retrieve genotype call and depths
            genotype = call['GT'] or '.'
            depth = int(call['DP'] or 0)
            allele_depths = call['AD'] or [0, 0, 0, 0]
        except:
            return record, "-", "-"
        # Ensure read counts match the alleles

        # Get all alleles for that record - start with the reference
        alleles = [record.REF] + record.ALT
        # Get the total read count, and test if we have enough coverage to make a call
        read_count = 100 #[ for index in range(0, )]
        if read_count < self.min_call_reads:
            return record, "-", "-"
        # If we have multiple alleles, filter out those alleles that have insufficient reads counts or read proportion
            
        # Update allele count and repeat the coverage test -- make this and the previous check a function
        
        if read_count < self.min_allele_reads and self.min_allele_proportion:
            return record, "-", "-"
        return record, alleles, depth

    def _format_record(self, record, genotypes, depth):
        '''
        Creates a row for output to a genotype file.
        '''
        genotypes_formatted = "".join(str(x) for x in genotypes).replace('None','')
        filter_value = record.FILTER[0] if len(record.FILTER) != 0 else 'PASS'
        return [record.CHROM, str(record.POS), str(genotypes_formatted), str(depth), str(filter_value)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf_list", "-i", help="List of VCF files to merge, mask SNPs in and update co-ordinates for.", required=True, type=str, nargs="+", default=[])
    parser.add_argument("--output_file_name", "-o", help="Base name (no extension) for the output genotypes file.", required=True)
    parser.add_argument("--chromKey_file", "-m", help="File containing SNPs to mask in merged VCF.", required=True)
    parser.add_argument("--MinCallReads", "-c", help=".", required=True)
    parser.add_argument("--MinAlleleReads", "-a", help=".", required=True)
    parser.add_argument("--MinAlleleProp", "-p", help=".", required=True)
    args = parser.parse_args()
    liftover_genotypes_write_file = GenotypeFileWriter(args.input_vcf_list, args.output_file_name, args.chromKey_file, args.MinCallReads, args.MinAlleleReads, args.MinAlleleProp)
    liftover_genotypes_write_file.write_genotype_file()
