import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy.metaphlan import BaseMetaphlanParser


class TestMetaphlan2Parser(TestCase):

    def setUp(self):
        input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.base_metaphlan_parser = BaseMetaphlanParser(input_path, analysis_type='marker_counts')

    def test_rows_differences_mismatch_counts(self):
        df1 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df1 = df1.set_index(['OTU ID'])

        df2 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 5.7],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df2 = df2.set_index(['OTU ID'])

        expected_df = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 0.0, 2.3]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['OTU ID'])

        observed_df = self.base_metaphlan_parser.rows_differences(df1, df2)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)

    def test_rows_differences_mismatch_counts_negative(self):
        df1 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 5.7],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df1 = df1.set_index(['OTU ID'])

        df2 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df2 = df2.set_index(['OTU ID'])

        expected_df = pd.DataFrame(
            [],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['OTU ID'])

        observed_df = self.base_metaphlan_parser.rows_differences(df1, df2)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_dtype=False)

    def test_rows_differences_missing_row(self):
        df1 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df1 = df1.set_index(['OTU ID'])

        df2 = pd.DataFrame(
            [
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df2 = df2.set_index(['OTU ID'])

        expected_df = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['OTU ID'])

        observed_df = self.base_metaphlan_parser.rows_differences(df1, df2)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)

    def test_rows_differences_missing_row_other_way(self):
        df1 = pd.DataFrame(
            [
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df1 = df1.set_index(['OTU ID'])

        df2 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df2 = df2.set_index(['OTU ID'])

        expected_df = pd.DataFrame(
            [],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['OTU ID'])

        observed_df = self.base_metaphlan_parser.rows_differences(df1, df2)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_dtype=False)

    def test_compare_difference_between_two_levels(self):
        df1 = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                 5, 3.6]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df1 = df1.set_index(['OTU ID'])

        df2 = pd.DataFrame(
            [
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|\
                  s__Actinobaculum_massiliense', 1.0, 2.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
                  s__Streptococcus_thermophilus', 1.7, 0.7],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
                  s__Streptococcus_salivarius', 3.3, 1.2]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        df2 = df2.set_index(['OTU ID'])

        expected_df = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                 0, 1.7]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )
        expected_df = expected_df.set_index(['OTU ID'])

        observed_df = self.base_metaphlan_parser.compare_difference_between_two_levels(df1, df2, 6)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)

    def test_remove_duplicates(self):

        tested_df = pd.DataFrame(
            [
                ['k__Bacteria', 10.5, 12.3],
                ['k__Bacteria|p__Actinobacteria', 1.0, 2.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria', 1.0, 2.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales', 1.0, 2.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae', 1.0, 2.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|\
g__Actinobaculum', 1.0, 2.0],
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|\
s__Actinobaculum_massiliense', 1.0, 2.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales', 9.5, 10.3],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae', 3.2, 8.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae', 6.3, 2.3],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                 6.3, 2.3],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_thermophilus', 1.7, 0.7],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_salivarius', 3.3, 1.2]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )

        expected_df = pd.DataFrame(
            [
                ['k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Actinomycetaceae|g__Actinobaculum|\
s__Actinobaculum_massiliense', 1.0, 2.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_thermophilus', 1.7, 0.7],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_salivarius', 3.3, 1.2],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus',
                 3.2, 8.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                 1.3, 0.4]

            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )

        observed_df = self.base_metaphlan_parser.remove_duplicates(tested_df)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)

    def test_remove_duplicates_rel_ab_addition_error_margin(self):
        input_path = os.path.join(os.path.dirname(__file__), 'input.tsv')
        self.base_metaphlan_parser = BaseMetaphlanParser(input_path, analysis_type='rel_ab')

        tested_df = pd.DataFrame(
            [
                ['k__Bacteria', 100.0, 100.0],
                ['k__Bacteria|p__Firmicutes', 100.0, 100.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli', 100.0, 100.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales', 100.0, 100.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae', 100.0, 100.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                 100.0, 100.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_thermophilus', 35.95876, 50.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_salivarius', 64.04123, 32.8]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )

        expected_df = pd.DataFrame(
            [
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_thermophilus', 35.95876, 50.0],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|\
s__Streptococcus_salivarius', 64.04123, 32.8],
                ['k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus',
                 0.0, 17.2]
            ],
            columns=['OTU ID', 'SAMPLE_1', 'SAMPLE_2']
        )

        observed_df = self.base_metaphlan_parser.remove_duplicates(tested_df)
        pd.testing.assert_frame_equal(observed_df, expected_df, check_like=True)
