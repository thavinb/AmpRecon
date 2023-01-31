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
    def __init__(self, vcf_file_path_list, out_file_name, chromKey_file_path):
        self.vcf_list = vcf_file_path_list
        self.output_file_name = out_file_name
        self.chromKey_file = chromKey_file_path

    def write_genotype_file(self):
        '''
        Merges supplied VCF files into a single genotype .tsv file.
        Replaces contigs field lines with the supplied contigs and their lengths.
        Updates VCF record co-ordinates to match those found within ChromKey data frame.
        SNPs to be masked are not written to the resulting file.
        '''
        chromKey_data_frame = pd.read_csv(self.chromKey_file, sep="\t", header=0)
        vcf_records = self._get_all_vcf_records(self.vcf_list)
        with open(self.output_file_name, "w") as output_genotype_file:
           output_genotype_file.write("Chr\tLoc\t\tGen\tIsHet\tDepth\tAle1_Dep\tAle2_Dep\tFilt\n")
           # Update and write data lines
           lifted_over_records = [row for row in [self._lift_over_record_coordinates(record, chromKey_data_frame) for record in vcf_records] if row != None]
           # Filtering and row formatting
           genotype_file_lines = [self._filter_and_format_row(line, chromKey_data_frame) for line in lifted_over_records]
           # Write rows to output genotype file
           [output_genotype_file.write(str(line)) for line in genotype_file_lines]

    def _get_all_vcf_records(self, vcf_list):
        '''
        Retrieves all the records found in all the supplied gzip compressed VCF files
        '''
        file_lines = []
        for vcf_file in vcf_list:
            vcf_file_open = vcf.Reader(filename=vcf_file)
            file_lines.extend(vcf_file_open)
        return file_lines

    def _lift_over_record_coordinates(self, vcf_record, chromKey_data_frame):
        '''
        Updates the chromosome and locus of a supplied VCF record to match those in the ChromKey data frame.
        VCF Records at positions to be masked are dropped.
        '''
        try:
            row = chromKey_data_frame.loc[(chromKey_data_frame["Chrom_ID"] == vcf_record.CHROM) & (chromKey_data_frame["VarPos"] == int(vcf_record.POS))]
            if row.iloc[0]["Mask"] == 0:
                vcf_record.CHROM = row.iloc[0]["Chromosome"]
                vcf_record.POS = row.iloc[0]["Locus"]
                return vcf_record
        except IndexError:
            logging.error(f"Failed to update data line co-ordinates - {vcf_row[0]}:{vcf_row[1]}")
            exit(1)

    def _filter_and_format_row(self, vcf_record, chromKey_data_frame):
        '''
        Filters and updates the VCF row allele depths and genotypes. Formats the rows for output to a genotype file.
        '''
        chromosome = vcf_record.CHROM
        locus = vcf_record.POS
#------Code below has been copied so far
        try:
            call = vcf_record.genotype(sample)
            genotypes = [call['GT'] or '.', 0, int(call['DP'] or 0), call['AD'] or [0, 0, 0, 0]]
            if not isinstance(genotypes[3], list):
                genotypes[3] = [genotypes[3]]
            while len(genotypes[3]) < 4:
                genotypes[3].append(0)
        except:
            call = None
            genotypes = ['0/0', 0, 0, [0, 0, 0, 0]]

        if vcf_record.ALT == [None]:
            actual_row = chromKey_data_frame.loc[(chromKey_data_frame["Chromosome"] == vcf_record.CHROM) & (chromKey_data_frame["Locus"] == vcf_record.POS)]
            actual = actual_row.iloc[0]["RefAllele"]
            if call is None:
                actual = '-'

            genotypes_out = actual
            depth = genotypes[2]
            allele_1_depth = genotypes[3][0]
            allele_2_depth = 0
        else:
            alleles = [vcf_record.REF] + vcf_record.ALT
            ad = genotypes[3]
            if genotypes[0] == '0/0':
                genotypes_out = str(alleles[0])
                depth = genotypes[2]
                allele_1_depth = ad[0]
                allele_2_depth = ad[1]
            elif genotypes[0] == '1/1':
                genotypes_out = str(alleles[1])
                depth = genotypes[2]
                allele_1_depth = ad[0]
                allele_2_depth = ad[1]
            elif genotypes[0] == '2/2':
                genotypes_out = str(alleles[2])
                depth = genotypes[2]
                allele_1_depth = genotypes[2] - ad[2]
                allele_2_depth = ad[2]
            else:
                one = int(genotypes[0][0])
                two = int(genotypes[0][2])
                genotypes_out = str(alleles[one]) + str(alleles[two])
                depth = genotypes[2]
                allele_1_depth = ad[one]
                allele_2_depth = ad[two]
#------Code above has been copied so far
        heterozygosity = '1' if len(genotypes_out) == 2 else '-'
        filter_value = vcf_record.FILTER[0] if len(vcf_record.FILTER) != 0 else 'PASS'
        columns = [chromosome, str(locus), genotypes_out, heterozygosity, str(depth), str(allele_1_depth), str(allele_2_depth), str(filter_value)]
        formatted_line = "\t".join(columns) + "\n"
        return formatted_line

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf_list", "-i", help="List of VCF files to merge, mask SNPs in and update co-ordinates for.", required=True, type=str, nargs="+", default=[])
    parser.add_argument("--output_file_name", "-o", help="Base name (no extension) for the output genotypes file.", required=True)
    parser.add_argument("--chromKey_file", "-m", help="File containing SNPs to mask in merged VCF.", required=True)
    args = parser.parse_args()
    liftover_genotypes_write_file = GenotypeFileWriter(args.input_vcf_list, args.output_file_name, args.chromKey_file)
    liftover_genotypes_write_file.write_genotype_file()
