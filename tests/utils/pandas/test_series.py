from unittest import TestCase

import pandas as pd

from moonstone.utils.pandas.series import SeriesStatsBuilder, SeriesBinning


class TestSeriesStatsBuilder(TestCase):

    def test_build_stats_str(self):
        series = pd.Series(
            ['M', 'F', 'M'],
            name='test_str'
        )
        expected_dict = {
            'col_name': 'test_str',
            'col_type': 'object',
            'python_col_type': 'str',
            'n_values': 3,
            'n_uniq_values': 2,
            'values_repartition': {
                'M': 2,
                'F': 1
            }
        }
        stats_builder = SeriesStatsBuilder(series)
        tested_dict = stats_builder.build_stats()
        self.assertDictEqual(tested_dict, expected_dict)

    def test_build_stats_int(self):
        series = pd.Series(
            [12, 41, 15],
            name='test_int'
        )
        expected_dict = {
            'col_name': 'test_int',
            'col_type': 'int64',
            'python_col_type': 'int',
            'n_values': 3,
            'n_uniq_values': 3,
            'mean': 22.67,
            'values_repartition': {
                12: 1,
                41: 1,
                15: 1,
            }
        }
        stats_builder = SeriesStatsBuilder(series)
        tested_dict = stats_builder.build_stats()
        self.assertDictEqual(tested_dict, expected_dict)

    def test_build_stats_float(self):
        series = pd.Series(
            [13.3, 15.3, 19.1],
            name='test_float'
        )
        expected_dict = {
            'col_name': 'test_float',
            'col_type': 'float64',
            'python_col_type': 'float',
            'n_values': 3,
            'n_uniq_values': 3,
            'mean': 15.9,
            'values_repartition': {
                13.3: 1,
                15.3: 1,
                19.1: 1,
            }
        }
        stats_builder = SeriesStatsBuilder(series)
        tested_dict = stats_builder.build_stats()
        self.assertDictEqual(tested_dict, expected_dict)


class TestSeriesBinning(TestCase):

    def test_compute_heterogeneous_bins(self):
        tested_object = pd.Series(
            {
                'gene_1': 10.5,
                'gene_2': 5.9,
                'gene_3': 9,
            })
        tested_object.name = 'mean read count'
        tested_object_instance = SeriesBinning(tested_object)
        tested_object_instance.bins_values     # call compute_heterogeneous_bins()
        expected_object = [-0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20.0]
        self.assertListEqual(tested_object_instance.bins_values, expected_object)

    def test_compute_binned_data(self):
        series = pd.Series(
            [5, 8, 10], index=['i1', 'i2', 'i3']
        )
        expected_object = pd.Series(
            [1, 2], index=[']0, 5]', ']5, 10]']
        )
        tested_object_instance = SeriesBinning(series)
        tested_object_instance.bins_values = [0, 5, 10]
        tested_object = tested_object_instance.compute_binned_data()
        pd.testing.assert_series_equal(tested_object, expected_object)
        pd.testing.assert_series_equal(tested_object_instance.binned_data, expected_object)
