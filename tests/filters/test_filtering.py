from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.filters.filtering import (
    Filtering
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

    def test_compute_threshold_best_n_percent(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [9, 7, 8, 3],
                'specie_4': [3, 2, 1, 0],
                'specie_5': [8, 3, 0, 1],
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object.columns.name = 'sample'
        tested_object_instance = Filtering(tested_object)
        expected_object = float(2.9)
        self.assertEqual(tested_object_instance.compute_threshold_best_n_percent(),
                         expected_object)

    def test_keep_data(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [9, 7, 8, 3],
                'specie_4': [3, 2, 1, 0],
                'specie_5': [8, 3, 0, 1],
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object.columns.name = 'sample'
        tested_object_instance = Filtering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [9, 7, 8, 3],
                'specie_5': [8, 3, 0, 1],
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.keep_data(),
                                      expected_object)
