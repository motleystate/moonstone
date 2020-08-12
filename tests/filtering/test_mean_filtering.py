from unittest import TestCase

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
