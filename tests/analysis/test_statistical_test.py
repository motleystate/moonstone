from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.analysis.statistical_test import (
    statistical_test_groups_comparison,
    _compute_best_bins_values
)


class TestStatisticalTestFunction(TestCase):

    def setUp(self):
        self.test_df = pd.Series({
            'sample1': 0.72,
            'sample2': 1.98,
            'sample3': 0.00,
            'sample4': 0.00,
            'sample5': 4.99
        })

    def test_mann_whitney_u_groups(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )

        expected_df = pd.DataFrame(
            [
                [np.nan, 0.374259],
                [0.374259, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u')
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_groups_wrong_metadata(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1]
            ],
            columns=['sex', 'group'],
            index=['sample6', 'sample7'],
        )

        with self.assertRaises(RuntimeError) as cm:
            statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u')
        the_exception = cm.exception
        self.assertEqual(the_exception.__str__(), "All groups have been dropped: not enough observations by group.")

    def test_mann_whitney_u_groups_nonsym(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )

        expected_df = pd.DataFrame(
            [
                [np.nan, 0.374259],
                [np.nan, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )
        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u',
                                                    output='dataframe', sym=False)
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_groups_series_output(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )

        expected_df = pd.Series({(1, 2): 0.3742593192802244, (2, 1): 0.3742593192802244})
        expected_df.index.names = ['Group1', 'Group2']

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u',
                                                    output='series', sym=True)
        pd.testing.assert_series_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_groups_nonsym_series_output(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )

        expected_df = pd.Series({(1, 2): 0.3742593192802244})
        expected_df.index.names = ['Group1', 'Group2']

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u',
                                                    output='series', sym=False)
        pd.testing.assert_series_equal(matrix, expected_df, check_dtype=False)

    def test_ttest_independance_groups(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )

        expected_df = pd.DataFrame(
            [
                [np.nan, 0.483379],
                [0.483379, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='ttest_independence')
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)


class TestChi2Functions(TestCase):

    def setUp(self):
        self.test_df = pd.Series([7, 11, 13, 18, 14, 17, 9, 0, 8, 17, 8, 19, 8, 11, 1, 6, 16, 16, 15, 18, 6, 4, 3, 4,
                                  18, 16, 8, 2, 7, 14, 11, 9, 15, 10, 9, 9, 1, 7, 14, 18, 14, 12, 18, 14, 8, 4, 3, 7,
                                  2, 5, 6, 3, 19, 19, 8, 12, 7, 1, 8, 1, 11, 14, 17, 13, 10, 6, 11, 1, 5, 18, 1, 11, 1,
                                  16, 17, 19, 9, 0, 17, 7, 18, 0, 15, 7, 2, 12, 7, 9, 1, 6, 15, 9, 5, 19, 13, 17, 1, 17,
                                  3, 0])
        self.metadata_df = pd.Series([3, 2, 1, 1, 1, 2, 3, 2, 2, 2, 2, 3, 1, 3, 1, 3, 2, 3, 1, 1, 3, 1, 2, 1, 1, 1, 1,
                                      3, 3, 2, 2, 2, 2, 2, 3, 2, 3, 3, 1, 2, 1, 3, 1, 2, 1, 1, 1, 1, 3, 3, 2, 2, 3, 3,
                                      2, 2, 2, 1, 1, 1, 1, 3, 1, 1, 3, 2, 2, 3, 1, 2, 1, 2, 1, 2, 3, 2, 2, 3, 3, 1, 3,
                                      3, 3, 1, 3, 2, 2, 2, 1, 2, 3, 2, 1, 1, 3, 2, 3, 1, 1, 2])

    def test_compute_best_bins_values(self):
        list_of_series = [
            self.test_df[self.metadata_df[self.metadata_df == 1].index].dropna(),
            self.test_df[self.metadata_df[self.metadata_df == 2].index].dropna(),
            self.test_df[self.metadata_df[self.metadata_df == 3].index].dropna()
        ]
        bins_values = _compute_best_bins_values(list_of_series)
        expected_bins_values = [-0.001, 6.333333333333333, 12.666666666666666, 19.0]
        self.assertListEqual(bins_values, expected_bins_values)

    def test_chi2_contingency_groups(self):
        expected_df = pd.DataFrame(
            [
                [np.nan, 0.06720551273974977, 0.9310443064194389],
                [0.06720551273974977, np.nan, 0.152090],
                [0.9310443064194389, 0.152090, np.nan]
            ],
            columns=[1, 2, 3],
            index=[1, 2, 3]
        )

        matrix = statistical_test_groups_comparison(self.test_df, self.metadata_df, stat_test='chi2_contingency')
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)
