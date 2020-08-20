from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.filtering.basics_filtering import (
    NoCountsFiltering, NamesFiltering
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
