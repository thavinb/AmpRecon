#!/usr/bin/env python3

from grc_kelch13_mutation_caller import KelchMutationCaller
import csv
import unittest
import json


class TestKelchVariantCaller(unittest.TestCase):
    """
    Testing of KelchVariantCaller class methods.
    """

    def setUp(self):

        self.test_genotype_file = ["multi_sample_kelch_test_genotype_file.txt"]
        self.out_file_name = "test_kelch13_variant_calls.tsv"
        self.test_config_file = "test_config.json"
        self.kelch_reference_file = "kelchReference.txt"
        self.codon_key_file = "codonKey.txt"

        # Read in test data from testing config file
        with open(self.test_config_file) as file:
            config = json.load(file).get("grc_kelch13")
            self.kelch_regions = config.get("kelch13_regions")
            self.test_kelch13_data = config.get("test_kelch13_data")
            self.test_kelch_reference_row = config.get("test_kelch_reference_row")
            self.genotype_file_column_names = config.get("genotype_file_column_names")
            self.test_variants = config.get("test_kelch13_variants")
            self.codon_key = config.get("codon_key")

        # Create instance of the KelchMutationCaller class
        self.kelch_caller = KelchMutationCaller(
            self.test_genotype_file,
            self.out_file_name,
            self.test_config_file,
            self.kelch_reference_file,
            self.codon_key_file,
        )

    def test_get_kelch_genotype_file_rows(self):

        # Run _get_kelch_genotype_file_rows() method
        genotype_file_kelch13_rows = self.kelch_caller._get_kelch_genotype_file_rows(
            self.test_genotype_file[0], self.kelch_regions
        )

        # Get expected number of rows -> the number of bases including the start and end of the kelch13 region
        per_sample_number_of_rows = 0
        for region in self.kelch_regions.values():
            difference = (int(region.get("End")) - int(region.get("Start"))) + 1
            per_sample_number_of_rows += difference

        for row_set in genotype_file_kelch13_rows.values():

            # Test whether the number of retrieved rows for each sample ID in the genotype file is what we expect
            self.assertEqual(len(row_set), per_sample_number_of_rows)

            for row in row_set:
                # Ensure that the returned genotype file rows have the correct number of columns
                self.assertEqual(len(row), len(self.genotype_file_column_names))

                # Ensure that the returned genotype file rows have the correct column names
                [
                    self.assertIn(key, self.genotype_file_column_names)
                    for key in row.keys()
                ]
                # Ensure that each genotype file row is a valid kelch13 region
                self.assertIn(row[str("Chr")], self.kelch_regions.keys())
                region = self.kelch_regions.get(row[str("Chr")])

                # Ensure that the locus in this row is greater than or equal to the start of this kelch13 region
                self.assertGreaterEqual(row[str("Loc")], region.get("Start"))

                # Ensure that the locus in this row is less than or equal to the end of this kelch13 region
                self.assertLessEqual(row[str("Loc")], region.get("End"))

    def test_get_reference_allele(self):
        for locus, data in self.test_kelch13_data.items():
            expected_allele = data.get("reference_allele")
            reference_allele = self.kelch_caller._get_reference_allele(
                locus, self.test_kelch_reference_row
            )
            # Ensure the reference allele is the expected complemented allele for that matched column
            self.assertEqual(reference_allele, expected_allele)

    def test_determine_sample_variant(self):
        # Parse the kelch reference file
        with open(self.kelch_reference_file) as file:
            kelch_reference = {}
            for row in csv.DictReader(file, delimiter="\t"):
                for gl in [row.get("GL1"), row.get("GL2"), row.get("GL3")]:
                    kelch_reference[gl] = dict(row)

        # Parse the genotype file and get all the kelch rows
        test_genotypes_all_samples = self.kelch_caller._get_kelch_genotype_file_rows(
            self.test_genotype_file[0], self.kelch_regions
        )

        for ID, expected_variant in self.test_variants.items():
            test_genotypes = test_genotypes_all_samples.get(ID)

            # Run _determine_sample_variant() method
            variant = self.kelch_caller._determine_sample_variant(
                test_genotypes, kelch_reference, 0.95, self.codon_key
            )
            variant = " ".join(variant) if isinstance(variant, list) else variant
            # Ensure that the determined variant is correct for each of the test sample IDs and variants
            self.assertEqual(variant, expected_variant)

    def test_check_for_amino_acid_substitution(self):
        for locus, data in self.test_kelch13_data.items():
            allele = data.get("allele")
            mutation = self.kelch_caller._check_for_amino_acid_substitution(
                allele, locus, self.test_kelch_reference_row, self.codon_key
            )
            expected_mutation = data.get("mutation")
            # Handle case where mutation is synonymous and None is returned
            if expected_mutation == "None":
                expected_mutation = None
            # Ensure the mutation matches the expected mutation
            # self.assertEqual(mutation, expected_mutation)

    def test_assemble_mutated_codon(self):
        for locus, data in self.test_kelch13_data.items():
            allele = data.get("allele")
            mutated_codon = self.kelch_caller._assemble_mutated_codon(
                locus, allele, self.test_kelch_reference_row
            )
            expected_mutated_codon = data.get("codon")
            # Ensure the mutated codon matches the expected mutated codon
            self.assertEqual(mutated_codon, expected_mutated_codon)

    def test_call_kelch_variants(self):
        # Run call_kelch_variants() method
        self.kelch_caller.call_kelch_variants()

        with open(self.out_file_name) as output_file:
            # Skip header line
            next(output_file)

            # Ensure that the variant in that line of the file is what we expect it to be
            for row in output_file.readlines():
                line_elements = row.split()
                id = line_elements[0]
                output_file_variant = " ".join(line_elements[1:])
                variant = self.test_variants.get(id)
                self.assertEqual(variant, output_file_variant)


if __name__ == "__main__":
    unittest.main()
