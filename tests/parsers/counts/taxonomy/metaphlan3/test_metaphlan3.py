import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import Metaphlan3Parser


class TestMetaphlan2Parser(TestCase):

    def setUp(self):
        self.input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.meta2parser = Metaphlan3Parser(self.input_path, analysis_type='marker_counts')

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

    def test_to_dataframe_keep_NCBI_tax_col(self):
        """
        Test based on input.tsv file
        """
        meta2parser = Metaphlan3Parser(self.input_path, analysis_type='rel_ab', keep_NCBI_tax_col=True)
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum',
                 'Actinobaculum_massiliense', 1.0, 2.0, '461393'],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus',
                 'Lactobacillus (genus)', 3.2, 8.0, '1578'],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 'Streptococcus (genus)', 1.3, 0.4, '1301'],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 'Streptococcus_thermophilus', 1.7, 0.7, '1308'],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 'Streptococcus_salivarius', 3.3, 1.2, '1304']
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'SAMPLE_1', 'SAMPLE_2',
                     'NCBI_tax_id']
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        observed_df = meta2parser.dataframe
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)

    def test_to_dataframe_less_taxonomical_names(self):
        """
        Test based on input.tsv file
        """
        meta2parser = Metaphlan3Parser(self.input_path, analysis_type='rel_ab', keep_NCBI_tax_col=True)
        meta2parser.taxonomical_names = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']
        expected_df = pd.DataFrame(
            [
                ['Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales', 'Actinomycetaceae', 'Actinobaculum',
                 1.0, 2.0, '76833'],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus',
                 3.2, 8.0, '1578'],
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
                 6.3, 2.3, '1301'],
            ],
            columns=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'SAMPLE_1', 'SAMPLE_2',
                     'NCBI_tax_id']
        )
        expected_df = expected_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
        observed_df = meta2parser.dataframe
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)
