#!/usr/bin/env python3

from plasmepsin_caller import PlasmepsinVariantCaller
import unittest
import os
import json


class TestPlasmepsinVariantCaller(unittest.TestCase):
    """
    Testing of PlasmepsinVariantCaller class methods.
    """

    def setUp(self):

        self.test_genotype_file = ["test_genotype_file.txt"]
        self.out_file_name = "test_plasmepsin_variant_calls.tsv"
        self.loci_genotypes_variants = "test_config.json"

        # Read in test data from testing config file
        with open(self.loci_genotypes_variants) as file:
            config = json.load(file)
            self.plasmepsin_loci = config.get("plasmepsin_loci")
            self.test_genotypes = config.get("test_genotypes")
            self.test_variants = config.get("test_variants")
            self.genotype_file_column_names = config.get("genotype_file_column_names")

        # Create instance of the PlasmepsinVariantCaller class
        self.plasmepsin_caller = PlasmepsinVariantCaller(
            self.test_genotype_file, self.out_file_name, self.loci_genotypes_variants
        )

    def test_get_plasmepsin_rows(self):
        plasmepsin_positions = [locus.get("Position") for locus in self.plasmepsin_loci]
        # Run _get_plasmepsin_rows() method
        genotype_file_plasmepsin_rows = self.plasmepsin_caller._get_plasmepsin_rows(
            self.test_genotype_file[0], self.plasmepsin_loci
        )
        for row in genotype_file_plasmepsin_rows.values():
            # Ensure that the returned genotype file rows have the correct number of columns
            self.assertEqual(len(row), len(self.genotype_file_column_names))

            # Ensure that the returned genotype file rows have the correct column names
            [self.assertIn(key, self.genotype_file_column_names) for key in row.keys()]

            # Ensure that each genotype file row is a plasmepsin locus
            position = f"{row[str('Chr')]}:{row[str('Chr_Loc')]}"
            self.assertIn(position, plasmepsin_positions)

    def test_determine_sample_variant(self):
        # Ensure that the determined variant is correct for each of the test sample IDs and genotypes
        for id, expected_variant in self.test_variants.items():
            # Run _determine_sample_variant() method
            variant = self.plasmepsin_caller._determine_sample_variant(
                id, self.test_genotypes, self.plasmepsin_loci
            )
            self.assertEqual(variant, expected_variant)

    def test_call_plasmepsin_variants(self):
        # Run call_plasmepsin_variants() method
        self.plasmepsin_caller.call_plasmepsin_variants()

        with open(self.out_file_name) as output_file:
            # Skip header line
            next(output_file)

            # Ensure that the variant in that line of the file is what we expect it to be
            for row in output_file.readlines():
                id, output_file_variant = row.split()
                variant = self.test_variants.get(id)
                self.assertEqual(variant, output_file_variant)


if __name__ == "__main__":
    unittest.main()
