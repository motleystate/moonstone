from unittest import TestCase

import pandas as pd

from moonstone.analysis.columns_statistics import DataframeStatistics


class TestDataframeStatistics(TestCase):

    def setUp(self):
        self.dataframe = pd.DataFrame([
            [1, "yes", 1.5],
            [4, "yes", 4.5],
            [1, "maybe", 3.0],
        ], columns=['int', 'choice', 'float'], index=['s1', 's2', 's3'])

    def test_get_stats_headers(self):

        expected_list = [
            {
                'col_name': 'int',
                'col_type': 'int64',
                'python_col_type': 'int',
                'n_values': 3,
                'n_uniq_values': 2,
                'mean': 2.0,
                'values_repartition': {
                    1: 2,
                    4: 1
                }
            },
            {
                'col_name': 'choice',
                'col_type': 'object',
                'python_col_type': 'str',
                'n_values': 3,
                'n_uniq_values': 2,
                'values_repartition': {
                    'yes': 2,
                    'maybe': 1
                }
            },
            {
                'col_name': 'float',
                'col_type': 'float64',
                'python_col_type': 'float',
                'n_values': 3,
                'n_uniq_values': 3,
                'mean': 3.0,
                'values_repartition': {
                    1.5: 1,
                    4.5: 1,
                    3.0: 1
                }
            }
        ]
        df_stats = DataframeStatistics(self.dataframe)
        self.assertListEqual(df_stats.get_stats(), expected_list)
