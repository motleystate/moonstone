import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import SunbeamKraken2Parser


class TestSunbeamKraken2Parser(TestCase):

    def setUp(self):
        input_path = os.path.join(os.path.dirname(__file__), 'sunbeam_kraken2.tsv')
        self.sunbeamkraken2parser = SunbeamKraken2Parser(input_path)

    def test_fill_none(self):
        taxa_df = pd.DataFrame(
            [
                ['Bacteria', 'Bacteroidetes'],
                ['Bacteria', None]
            ],
            columns=['kingdom', 'phylum']
        )
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Bacteroidetes'],
                ['Bacteria', 'Bacteria (kingdom)']
            ],
            columns=['kingdom', 'phylum']
        )
        tested_df = self.sunbeamkraken2parser._fill_none(taxa_df)
        pd.testing.assert_frame_equal(tested_df, expected_df)

    def test_taxa_fill_none(self):
        df = pd.DataFrame(
            [
                [2, 5.5, 6.5, 'k__Bacteria; p__; c__; o__; f__; g__; s__'],
                [109790, 9.5, 7.3, 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__jensenii'],  # noqa
            ],
            columns=[self.sunbeamkraken2parser.new_otu_id_name, 'sample_1', 'sample_2', self.sunbeamkraken2parser.taxa_column]  # noqa
        )
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 2, 5.5, 6.5],  # noqa
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus', 'jensenii', 109790, 9.5, 7.3]  # noqa
            ],
            columns=[
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
                self.sunbeamkraken2parser.new_otu_id_name, 'sample_1', 'sample_2'
            ]
        )
        tested_df = self.sunbeamkraken2parser.split_taxa_fill_none(df)
        print(tested_df)
        pd.testing.assert_frame_equal(tested_df, expected_df)

    def test_to_dataframe(self):
        """
        Test based on input.tsv file
        """
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 'Bacteria (kingdom)', 2, 5.5, 6.0],  # noqa
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillales (order)', 'Lactobacillales (order)', 'Lactobacillales (order)', 186826, 4.3, 2.1],  # noqa
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus', 'jensenii', 109790, 1.0, 2.0]  # noqa
            ],
            columns=[
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
                self.sunbeamkraken2parser.new_otu_id_name, 'SAMPLE_1', 'SAMPLE_2'
            ]
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        pd.testing.assert_frame_equal(self.sunbeamkraken2parser.dataframe, expected_df)
