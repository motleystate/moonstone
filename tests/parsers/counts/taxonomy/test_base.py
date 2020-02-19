from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy.base import TaxonomyCountsBaseParser


class TestTaxonomyCountsBaseParser(TestCase):

    def setUp(self):
        self.parser = TaxonomyCountsBaseParser('unexisting_file')

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
        tested_df = self.parser._fill_none(taxa_df)
        pd.testing.assert_frame_equal(tested_df, expected_df)

    def test_taxa_fill_none(self):
        df = pd.DataFrame(
            [
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae', 5.5, 6.5],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales', 9.5, 7.3],
            ],
            columns=[self.parser.taxa_column, 'sample_1', 'sample_2']
        )
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 5.5, 6.5],
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetales (order)', 9.5, 7.3]
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'sample_1', 'sample_2']
        )
        tested_df = self.parser.split_taxa_fill_none(df, sep='|')
        pd.testing.assert_frame_equal(tested_df, expected_df)
