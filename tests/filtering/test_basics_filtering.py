from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.filtering.basics_filtering import (
    Filtering, NoCountsFiltering, NamesFiltering
)


class TestFiltering(TestCase):

    def test_remove_data_without_metadata(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object.columns.name = 'sample'
        tested_object_metadata = pd.DataFrame.from_dict(
            {
                'Sex': ['M', 'F', 'M'],
                'AGE': [25, 65, 49]
            },
            orient='index', columns=['1', '3', '4'])
        tested_object_metadata.columns.name = 'sample'
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 1, 0],
                'specie_2': [25, 3, 9]
            },
            orient='index', columns=['1', '3', '4'])
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.remove_data_without_metadata(tested_object_metadata),
                                      expected_object)

    def test_remove_rows_wihtout_relevant_info(self):
        tested_object = pd.DataFrame.from_dict(
            {
                ('genus1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        tested_object.columns.name = 'sample'
        row_to_remove = [np.nan]
        level_to_consider = 0
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                ('genus1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9}
            },
            orient='index')
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.remove_rows_without_relevant_info(row_to_remove,
                                      level_to_consider), expected_object)

    def test_remove_rows_without_relevant_info_custom_value(self):
        tested_object = pd.DataFrame.from_dict(
            {
                ('genus1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'custom_value'): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        tested_object.columns.name = 'sample'
        row_to_remove = ['custom_value']
        # Test to filter out when we see custom_value in level 1
        level_to_consider = 1
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                ('genus1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9}
            },
            orient='index')
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.remove_rows_without_relevant_info(row_to_remove,
                                      level_to_consider), expected_object)
        # Verify that filtering based on level 0 raise an error
        level_to_consider = 0
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                ('genus1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'custom_value'): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        expected_object.columns.name = 'sample'
        with self.assertRaises(KeyError):
            pd.testing.assert_frame_equal(tested_object_instance.remove_rows_without_relevant_info(row_to_remove,
                                          level_to_consider), expected_object)

    def test_selecting_certain_rows(self):
        tested_object = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        tested_object.columns.name = 'sample'
        tested_object.index.set_names(['genus', 'specie'], inplace=True)
        desired_row_series = ['specie_1', 'specie_2']
        level_to_consider = 1
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9}
            },
            orient='index')
        expected_object.index.set_names(['genus', 'specie'], inplace=True)
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.selecting_rows(desired_row_series,
                                      level_to_consider), expected_object)

    def test_deleting_only_zeros_rows(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 0, 0, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object.columns.name = 'sample'
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.deleting_only_zeros_rows(tested_object),
                                      expected_object)


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
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=True)
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
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=1, keep=False)
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
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=0, keep=True)
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
        tested_filtering = NamesFiltering(test_df, selected_rows, axis=0, keep=False)
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
