from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.filtering.mean_filtering import (
    MeanFiltering
)


class TestMeanFiltering(TestCase):

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
        tested_object_instance = MeanFiltering(tested_object)
        expected_object = float(2.9)
        self.assertEqual(tested_object_instance.compute_threshold_best_n_percent(),
                         expected_object)

    def test_filter_no_threshold_given(self):
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
        tested_object_instance = MeanFiltering(tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [9, 7, 8, 3],
                'specie_5': [8, 3, 0, 1],
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.filter(),
                                      expected_object)

    def test_filter_threshold_given(self):
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
        tested_object_instance = MeanFiltering(tested_object, threshold=3.6)
        expected_object = pd.DataFrame.from_dict(
            {
                'specie_2': [25, 6, 3, 9],
                'specie_3': [9, 7, 8, 3],
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.filter(),
                                      expected_object)

    def test_filter_threshold_zero(self):
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
        tested_object_instance = MeanFiltering(tested_object, threshold=0)
        expected_object = tested_object.copy()
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object_instance.filter(),
                                      expected_object)

    def test_generate_report_data_no_threshold_given(self):
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
        tested_object_instance = MeanFiltering(tested_object)
        tested_object_instance.filter()
        expected_object = {
            'title': 'Filtering by mean',
            'data': {
                'threshold': 2.9,
                'n_items_removed': 1,
                'n_reads_removed': 6,
                'percentage_to_keep': 90
                }
            }
        self.assertDictEqual(tested_object_instance.generate_report_data(), expected_object)

    def test_generate_report_data_threshold_given(self):
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
        tested_object_instance = MeanFiltering(tested_object, 3.6)
        tested_object_instance.filter()
        # tested_object_instance.visualize("file_testing_visualization_filtering_by_mean.html")
        expected_object = {
            'title': 'Filtering by mean',
            'data': {
                'threshold': 3.6,
                'n_items_removed': 3,
                'n_reads_removed': 30
                }
            }
        self.assertDictEqual(tested_object_instance.generate_report_data(), expected_object)

    def test_mean_equal_to_threshold_kept(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [2, 2],
                'specie_2': [1, 2],
            },
            orient='index', columns=['1', '2'])
        filtered_df = MeanFiltering(tested_object, threshold=2).filtered_df
        self.assertEqual(list(filtered_df.index), ['specie_1'])

    def test_input_not_modified(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [0, 0, 1, 0],
            },
            orient='index', columns=['1', '2', '3', '4'])
        original = tested_object.copy()
        filtered_df = MeanFiltering(tested_object, threshold=1).filtered_df
        self.assertEqual(list(filtered_df.index), ['specie_1'])
        pd.testing.assert_frame_equal(tested_object, original)

    def test_row_with_no_values_removed(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [5.0, 5.0],
                'specie_2': [np.nan, np.nan],
            },
            orient='index', columns=['1', '2'])
        with self.assertLogs('moonstone.filtering.mean_filtering', level='WARNING') as logs:
            filtered_df = MeanFiltering(tested_object, threshold=1).filtered_df
        self.assertEqual(list(filtered_df.index), ['specie_1'])
        self.assertIn('specie_2', logs.output[0])

    def test_text_counts_raise_clear_error(self):
        tested_object = pd.DataFrame({'1': ['10', '12,5'], '2': [1.0, 2.0]}, index=['specie_1', 'specie_2'])
        with self.assertRaisesRegex(ValueError, 'Counts must be numbers.*1'):
            MeanFiltering(tested_object, threshold=1)

    def test_negative_counts_warn(self):
        tested_object = pd.DataFrame({'1': [-1.0, 2.0], '2': [1.0, 2.0]}, index=['specie_1', 'specie_2'])
        with self.assertLogs('moonstone.filtering.base', level='WARNING') as logs:
            MeanFiltering(tested_object, threshold=1)
        self.assertIn('1 negative count', logs.output[0])
