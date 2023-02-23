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
        out_file_name = "test_plasmepsin_variant_calls.tsv"
        self.loci_genotypes_variants = "test_config.json"

        with open(self.loci_genotypes_variants) as file:
            config = json.load(file)
            self.plasmepsin_loci = config.get("plasmepsin_loci")
            self.test_genotypes = config.get("test_genotypes")
            self.genotype_file_column_names = config.get("genotype_file_column_names")

        self.plasmepsin_caller = PlasmepsinVariantCaller(
            self.test_genotype_file, out_file_name, self.loci_genotypes_variants
        )

    def test_get_plasmepsin_rows(self):
        plasmepsin_positions = [locus.get("Position") for locus in self.plasmepsin_loci]
        genotype_file_plasmepsin_rows = self.plasmepsin_caller._get_plasmepsin_rows(
            self.test_genotype_file[0], self.plasmepsin_loci
        )
        # Ensure that the returned genotype file rows have the correct columns
        for row in genotype_file_plasmepsin_rows.values():
            self.assertEqual(len(row), len(self.genotype_file_column_names))
            [self.assertIn(key, self.genotype_file_column_names) for key in row.keys()]

            # Ensure that position is a plasmepsin one
            position = f"{row[str('Chr')]}:{row[str('Chr_Loc')]}"
            self.assertIn(position, plasmepsin_positions)

    def test_call_variant(self):
        # Ensure that the determined variant is correct for each of the test sample IDs and genotypes
        sample_ids = ["ID_1", "ID_2", "ID_3", "ID_4"]
        test_variants = ["Copy", "WT", "-", "Copy"]

        for sample_id, test_variant in zip(sample_ids, test_variants):
            variant = self.plasmepsin_caller._call_variant(
                sample_id, self.test_genotypes, self.plasmepsin_loci
            )
            self.assertEqual(variant, test_variant)


if __name__ == "__main__":
    unittest.main()
