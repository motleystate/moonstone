from unittest import TestCase

import pandas as pd

from moonstone.utils.pandas.series import SeriesStatsBuilder


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
