import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import (
    Qiime2Parser
)


class TestQiime2Parser(TestCase):

    def setUp(self):
        input_path = os.path.join(os.path.dirname(__file__), 'raw_qiime2.csv')
        self.qiime2parser = Qiime2Parser(input_path, sep=',')

    def test_to_dataframe(self):
        """
        Test based on input.tsv file
        """
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Bacteroidetes', 'Bacteroidia', 'Bacteroidales', 'Tannerellaceae', 'Macellibacteroides', 0, 0, 0, 1],  # noqa
                ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Betaproteobacteriales', 'Chitinibacteraceae', 'Chitinibacteraceae (family)', 0, 2, 0, 0],  # noqa
                ['Bacteria', 'Actinobacteria', 'Acidimicrobiia', 'Microtrichales', 'Microtrichaceae', 'Microtrichaceae (family)', 0, 0, 3, 0],  # noqa
                ['Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Caulobacterales', 'Hyphomonadaceae', 'Hyphomonadaceae (family)', 5, 0, 0, 0],  # noqa
            ],
            columns=[
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                'Sample 1', 'Sample 2', 'Sample 3', 'Sample 4'
            ]
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
        pd.testing.assert_frame_equal(self.qiime2parser.dataframe, expected_df)
