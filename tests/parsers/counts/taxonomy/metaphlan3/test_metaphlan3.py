import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import Metaphlan3Parser


class TestMetaphlan2Parser(TestCase):

    def setUp(self):
        input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.meta2parser = Metaphlan3Parser(input_path)

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

