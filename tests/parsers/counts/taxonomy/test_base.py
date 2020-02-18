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
