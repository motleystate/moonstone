from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.filtering.basics_filtering import (
    NoCountsFiltering, NamesFiltering, NaNPercentageFiltering,
    NumberOfDifferentValuesFiltering
)


class TestNoCountsFiltering(TestCase):

    def test_int_df_rows_filtering(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 0, 0, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        tested_filtering = NoCountsFiltering(test_df, axis=1)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_float_df_rows_filtering(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3.0, 2.0, 1.0, 0.0],
                'specie_2': [25.0, 6.0, 3.0, 9.0],
                'specie_3': [0.0, 0.0, 0.0, 0.0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        tested_filtering = NoCountsFiltering(test_df, axis=1)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3.0, 2.0, 1.0, 0.0],
                'specie_2': [25.0, 6.0, 3.0, 9.0],
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_int_df_columns_filtering(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 0, 1, 0],
                'specie_2': [25, 0, 3, 9],
                'specie_3': [0, 0, 0, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        tested_filtering = NoCountsFiltering(test_df, axis=0)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 1, 0],
                'specie_2': [25, 3, 9],
                'specie_3': [0, 0, 0]
            },
            orient='index', columns=['1', '3', '4'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)


class TestNamesFiltering(TestCase):

    def test_selecting_rows(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        selected_rows = ['specie_1', 'specie_2']
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=0, keep=True)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_removing_rows(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        selected_rows = ['specie_1', 'specie_2']
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=0, keep=False)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_selecting_columns(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        selected_rows = ['1', '3']
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=True)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 1],
                'specie_2': [25, 3],
                'specie_3': [0, 5]
            },
            orient='index', columns=['1', '3'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_removing_columns(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        selected_rows = ['1', '3']
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=False)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [2, 0],
                'specie_2': [6, 9],
                'specie_3': [7, 0]
            },
            orient='index', columns=['2', '4'])
        expected_df.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_selecting_rows_multiindex(self):
        test_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        test_df.columns.name = 'sample'
        test_df.index.set_names(['genus', 'species'], inplace=True)
        selected_rows = ['specie_1', 'specie_2']
        with self.assertRaises(TypeError):
            tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=True)  # noqa

    def test_selecting_rows_without_names_present(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        selected_rows = ['specie_1', 'specie_5']
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        with self.assertLogs('moonstone.filtering.basics_filtering', level='INFO') as log:
            tested_filtering = NamesFiltering(test_df, selected_rows, axis=0, keep=True)
            pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)
            self.assertEqual(len(log.output), 1)
            self.assertIn("INFO:moonstone.filtering.basics_filtering:['specie_5']: \
1 name(s) not found in the dataframe.", log.output)

    def test_selecting_columns_without_names_present(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        selected_rows = ['1', '5', '6']
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=True)
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3],
                'specie_2': [25],
                'specie_3': [0]
            },
            orient='index', columns=['1'])
        expected_df.columns.name = 'sample'
        with self.assertLogs('moonstone.filtering.basics_filtering', level='INFO') as log:
            tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=True)
            pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)
            self.assertEqual(len(log.output), 1)
            self.assertIn("INFO:moonstone.filtering.basics_filtering:['5', '6']: \
2 name(s) not found in the dataframe.", log.output)


class TestByPercentageNaNFiltering(TestCase):

    def setUp(self):
        self.test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [np.nan, 2, 1, 0],
                'specie_2': [np.nan, 6, np.nan, np.nan],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        self.test_df.columns.name = 'sample'

    def test_filter_rows(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [np.nan, 2, 1.0, 0.0],         # transform into float during filtering -> no idea why
                'specie_3': [0.0, 7, 5.0, 0.0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=50, axis=0)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_filter_columns(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [2, 1, 0],
                'specie_2': [6, np.nan, np.nan],
                'specie_3': [7, 5, 0]
            },
            orient='index', columns=['2', '3', '4'])
        expected_df.columns.name = 'sample'
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=50, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_filter_columns_70_nan(self):
        expected_df = self.test_df
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=70, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_filter_rows_25_nan(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_1': [np.nan, 2, 1, 0],
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=25, axis=0)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)

    def test_filter_rows_0_nan(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'specie_3': [0, 7, 5, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_df.columns.name = 'sample'
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=0, axis=0)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)

    def test_filter_rows_100_nan(self):
        expected_df = self.test_df
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=100, axis=0)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)

    def test_fiter_fraction(self):
        expected_df = self.test_df
        tested_filtering = NaNPercentageFiltering(self.test_df, percentage=200/3, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)


class TestNumberOfDifferentValuesFiltering(TestCase):

    def setUp(self):
        self.test_df = pd.DataFrame.from_dict(
            {
                'sample_1': ['F', 35, 'n', 60],
                'sample_2': ['F', 19, np.nan, np.nan],
                'sample_3': ['F', 27, 'n', 55]
            },
            orient='index', columns=['sex', 'age', 'smoker', 'HDL cholesterol'])
        self.test_df.columns.name = 'metadata'

    def test_filter_not_counting_nan_columns(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'sample_1': [35, 60],
                'sample_2': [19, np.nan],
                'sample_3': [27, 55]
            },
            orient='index', columns=['age', 'HDL cholesterol'])
        expected_df.columns.name = 'metadata'
        tested_filtering = NumberOfDifferentValuesFiltering(self.test_df, min=2, na=False, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_filter_counting_nan_columns(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'sample_1': [35, 'n', 60],
                'sample_2': [19, np.nan, np.nan],
                'sample_3': [27, 'n', 55]
            },
            orient='index', columns=['age', 'smoker', 'HDL cholesterol'])
        expected_df.columns.name = 'metadata'
        tested_filtering = NumberOfDifferentValuesFiltering(self.test_df, min=2, na=True, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)

    def test_filter_rows(self):
        test_df_rows = self.test_df.transpose()
        expected_df = pd.DataFrame.from_dict(
            {
                'age': [35, 19, 27],
                'HDL cholesterol': [60, np.nan, 55],
            },
            orient='index', columns=['sample_1', 'sample_2', 'sample_3'])
        expected_df.index.name = 'metadata'
        tested_filtering = NumberOfDifferentValuesFiltering(test_df_rows, min=2, na=False, axis=0)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)

    def test_filter_max(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'sample_1': ['F', 'n', 60],
                'sample_2': ['F', np.nan, np.nan],
                'sample_3': ['F', 'n', 55]
            },
            orient='index', columns=['sex', 'smoker', 'HDL cholesterol'])
        expected_df.columns.name = 'metadata'
        tested_filtering = NumberOfDifferentValuesFiltering(self.test_df, max=2, na=False, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df, check_dtype=False)

    def test_filter_min_max(self):
        expected_df = pd.DataFrame.from_dict(
            {
                'sample_1': [60],
                'sample_2': [np.nan],
                'sample_3': [55]
            },
            orient='index', columns=['HDL cholesterol'])
        expected_df.columns.name = 'metadata'
        tested_filtering = NumberOfDifferentValuesFiltering(self.test_df, min=2, max=2, na=False, axis=1)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)
