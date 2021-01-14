import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import Metaphlan2Parser


class TestMetaphlan2Parser(TestCase):

    def test_to_dataframe(self):
        input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.meta2parser = Metaphlan2Parser(input_path, analysis_type='marker_counts')
        """
        Test based on input.tsv file
        """
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinomycetaceae (family)', 'Actinomycetaceae (family)', 1.2, 3.9],  # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum (genus)', 3.3, 0.1],   # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum_massiliense', 1.0, 2.0]   # noqa
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        pd.testing.assert_frame_equal(self.meta2parser.dataframe, expected_df, check_like=True)

    def test_to_dataframe_with_t__(self):
        """
        Test based on input_with_t__.tsv file, meaning it has annotation until the t__ level
        """
        input_path = os.path.join(os.path.dirname(__file__), 'input_with_t__.tsv')
        self.meta2parser = Metaphlan2Parser(input_path, analysis_type='marker_counts')
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinomycetaceae (family)', 'Actinomycetaceae (family)', 'Actinomycetaceae (family)', 1.2, 3.9],  # noqa
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum (genus)', 'Actinobaculum (genus)', 3.3, 0.1],   # noqa
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
        pd.testing.assert_frame_equal(self.meta2parser.dataframe, expected_df, check_like=True)
