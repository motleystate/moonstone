import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import Metaphlan2Parser


class TestMetaphlan2Parser(TestCase):

    def setUp(self):
        input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.meta2parser = Metaphlan2Parser(input_path)

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
        tested_df = self.meta2parser._fill_none(taxa_df)
        pd.testing.assert_frame_equal(tested_df, expected_df)

    def test_taxa_fill_none(self):
        df = pd.DataFrame(
            [
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae', 5.5, 6.5],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales', 9.5, 7.3],
            ],
            columns=[self.meta2parser.taxa_column, 'sample_1', 'sample_2']
        )
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 5.5, 6.5],
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetales (order)', 9.5, 7.3]
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'sample_1', 'sample_2']
        )
        tested_df = self.meta2parser.split_taxa_fill_none(df)
        pd.testing.assert_frame_equal(tested_df, expected_df)

    def test_to_dataframe(self):
        """
        Test based on input.tsv file
        """
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinomycetaceae (family)', 'Actinomycetaceae (family)', 5.5, 6.0],  # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum (genus)', 4.3, 2.1],   # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum_massiliense', 1.0, 2.0]   # noqa
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        pd.testing.assert_frame_equal(self.meta2parser.dataframe, expected_df)

    def test_to_dataframe_with_t__(self):
        """
        Test based on input_with_t__.tsv file, meaning it has annotation until the t__ level
        """
        input_path = os.path.join(os.path.dirname(__file__), 'input_with_t__.tsv')
        self.meta2parser = Metaphlan2Parser(input_path)
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinomycetaceae (family)', 'Actinomycetaceae (family)', 'Actinomycetaceae (family)', 5.5, 6.0],  # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum (genus)', 'Actinobaculum (genus)', 4.3, 2.1],   # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum_massiliense', 'GCF_000315465', 1.0, 2.0]   # noqa
            ],
            columns=[
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'sTrain', 'SAMPLE_1', 'SAMPLE_2'
            ]
        )
        expected_df = expected_df.set_index(
            [
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'sTrain'
            ]
        )
        pd.testing.assert_frame_equal(self.meta2parser.dataframe, expected_df)
