#!/usr/bin/env python3

import argparse
import logging
import vcf
import csv

class GenotypeFileWriter:
    '''
    Merges a list of VCF files into a single genotype .tsv file.
    Drops records at specified loci, updates the co-ordinates of the rows and filters out alleles with low coverage.
    '''
    def __init__(self, vcf_path_list, out_file_name, sample_id, chromKey_file_path, chromosome_column, locus_column, min_total_depth, het_min_allele_depth, het_min_allele_proportion):
        self.vcf_list = vcf_path_list
        self.output_file_name = out_file_name
        self.sample_id = sample_id
        self.chromKey_file = chromKey_file_path
        self.chromosome_column = chromosome_column
        self.locus_column = locus_column
        self.min_total_depth = float(min_total_depth)
        self.het_min_allele_depth = float(het_min_allele_depth)
        self.het_min_allele_proportion = float(het_min_allele_proportion)

    def write_genotype_file(self):
        '''
        Merges supplied VCF files into a single genotype .tsv file.
        Updates VCF record co-ordinates to match those found within ChromKey dictionary.
        SNPs to be masked are not written to the resulting file.
        '''
        # Read ChromKey file rows as dictionaries
        chromKey_dict = {}
        with open(self.chromKey_file, newline='') as chromKey_file:
            for row in csv.DictReader(chromKey_file, delimiter='\t'):
                key = f"{row[str(self.chromosome_column)]}:{row[str(self.locus_column)]}"
                chromKey_dict[key] = dict(row)

        # Prepare output genotype file
        output_genotype_file = open(self.output_file_name, "w")
        output_genotype_file.write("Sample_ID\tAmplicon\tAmplicon_Pos\tChr\tChr_Loc\tGen\tDepth\tFilt\n")

        genotype_file_rows = []

        # Iterate over the supplied VCF files and their records
        for vcf_file in self.vcf_list:
            vcf_reader = vcf.Reader(filename=vcf_file)
            for record in vcf_reader:

                # Match - Get row associated with record from chromKey file
                chromKey_row = self._match_chromKey_row(record, chromKey_dict)

                # Mask - Ignore this record if its supposed to be masked
                if chromKey_row.get("Mask") == "1":
                    continue

                # Lift over - Update record co-ordinates
                record, amplicon = self._lift_over_record_coordinates(record, chromKey_row)

                # Filter - Remove low coverage genotypes and adjust depths
                genotype, depth = self._filter_genotypes(record, vcf_reader)

                # Format record
                sample_id = self.sample_id or vcf_reader.samples[0]
                genotype_file_rows.append(self._format_record(sample_id, amplicon, record, genotype, depth))

        # Format each chromKey row and add any VCF record data at that position
        genotypes_dict =  {f"{row[1]}:{row[2]}" : row for row in genotype_file_rows}
        formatted_rows = [row for row in [self._format_chromKey_row(ck_row, genotypes_dict) for ck_row in chromKey_dict.items()] if row != None]

        # Write formatted chromKey / VCF record rows to genotype file
        output_genotype_file.writelines(formatted_rows)
        output_genotype_file.close()

    def _match_chromKey_row(self, record, chromKey_dict):
        '''
        Matches a VCF record to its associated dictionary in the chromKey dictionary.
        Matches by supplied columns (default is the amplicon and amplicon position)
        '''
        try:
            matching_chromKey_row = chromKey_dict.get(f"{record.CHROM}:{record.POS}")
            return matching_chromKey_row
        except IndexError:
            logging.error(f"Failed retrieve matching row from chromKey file for record co-ordinates: {record.CHROM}:{record.POS}")
            exit(1)

    def _lift_over_record_coordinates(self, record, chromKey_row):
        '''
        Updates the chromosome and locus of a supplied VCF record to match those in an associated dictionary.
        Also returns amplicon and amplicon position.
        '''
        amplicon = [chromKey_row.get("Chrom_ID"), chromKey_row.get("VarPos")] 
        record.CHROM = chromKey_row.get("Chromosome")
        record.POS = chromKey_row.get("Locus")
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
        Return alleles and depths that exceed depth and proportion thresholds.
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

    def _format_record(self, sample_id, amplicon, record, genotypes, allele_depths):
        '''
        Creates a row for output to a genotype file.
        '''
        genotypes_formatted = ",".join(str(genotype) for genotype in genotypes)
        allele_depths_formatted = ",".join(str(depth) for depth in allele_depths) if isinstance(allele_depths, list) else allele_depths
        filter_value = ";".join(record.FILTER) if record.FILTER is not None and len(record.FILTER) > 0 else "PASS"
        return [str(sample_id), amplicon[0], amplicon[1], record.CHROM, str(record.POS), str(genotypes_formatted), str(allele_depths_formatted), str(filter_value)]

    def _format_chromKey_row(self, ck_row, genotypes_dict):
        '''
        Format a chromKey row and adds any VCF record data at that position
        '''
        # Ensure position is not to be masked
        if ck_row[1].get("Mask") == '1':
            return None

        # If a genotype row taken from a VCF matches then use its values, otherwise use "-"
        genotypes_row = genotypes_dict.get(ck_row[0])
        if genotypes_row != None:
            genotype = genotypes_row[5]
            depth = genotypes_row[6]
            filter = genotypes_row[7]
            sample_id = genotypes_row[0]
        else:
            genotype = "-"
            depth = "-"
            filter = "-"
            sample_id = "-"

        # Retrieve location data from this chromKey row
        amplicon = ck_row[1].get("Chrom_ID")
        amplicon_pos = ck_row[1].get("VarPos")
        chr = ck_row[1].get("Chromosome")
        chr_loc = ck_row[1].get("Locus")

        # Format row for the output file
        row = [sample_id, amplicon, str(amplicon_pos), chr, str(chr_loc), genotype, depth, filter]
        return "\t".join(row)+"\n"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf_list", "-i", help="List of VCF files to merge, mask SNPs in and update co-ordinates for.", required=True, type=str, nargs="+", default=[])
    parser.add_argument("--output_file_name", "-o", help="Base name (no extension) for the output genotypes file.", type=str, required=True)
    parser.add_argument("--sample_id", "-s", help="Sample ID of the VCF files - output to genotype file.")
    parser.add_argument("--chromKey_file", "-m", help="File containing SNPs to mask in merged VCF.", type=str, required=True)
    parser.add_argument("--chromosome_column_name", "-c", help="Name of the column in the chromKey file to try to match chromosome to.", type=str, default="Chrom_ID")
    parser.add_argument("--locus_column_name", "-l", help="Name of the column in the chromKey file to try to match position to.", type=str, default="VarPos")
    parser.add_argument("--min_total_depth", "-d", help="The number of reads at a record must exceed this value.", type=float, default=10)
    parser.add_argument("--het_min_allele_depth", "-a", help="The number of reads for a particular allele at a heterozygous record must exceed this value.", type=float, default=5)
    parser.add_argument("--het_min_allele_proportion", "-p", help="The proportion of reads for particular allele at a heterozygous record must exceed this value.", type=float, default=0.1)
    args = parser.parse_args()
    liftover_genotypes_write_file = GenotypeFileWriter(args.input_vcf_list, args.output_file_name, args.sample_id, args.chromKey_file, args.chromosome_column_name, args.locus_column_name, args.min_total_depth, args.het_min_allele_depth, args.het_min_allele_proportion)
    liftover_genotypes_write_file.write_genotype_file()
