import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import Metaphlan3Parser


class TestMetaphlan2Parser(TestCase):

    def setUp(self):
        input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.meta2parser = Metaphlan3Parser(input_path, analysis_type='marker_counts')

    def test_to_dataframe(self):
        """
        Test based on input.tsv file
        """
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum',
                 'Actinobaculum_massiliense', 1.0, 2.0],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus',
                 'Lactobacillus (genus)', 3.2, 8.0],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 'Streptococcus (genus)', 1.3, 0.4],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 'Streptococcus_thermophilus', 1.7, 0.7],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 'Streptococcus_salivarius', 3.3, 1.2]
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        pd.testing.assert_frame_equal(self.meta2parser.dataframe, expected_df, check_like=True)
