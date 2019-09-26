import os
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from metautils.parsers.metadata import MetadataParser


class TestMetadataParser(TestCase):

    def setUp(self):
        self.test_file = os.path.join(os.path.dirname(__file__), "metadata.tsv")
        self.parser = MetadataParser(self.test_file)

    def test_parse_file(self):
        expected_df = pd.DataFrame(
            {
                'col_1': ['s1', 's2', 's3'],
                'col_2': [13.3, 15.3, 19.1],
                'col_3': ['M', 'F', 'M']
            }
        )
        assert_frame_equal(self.parser.dataframe, expected_df)

    def test_get_stats(self):
        expected_list = [
            {
                'col_name': 'col_1',
                'col_type': 'object',
                'python_col_type': 'str',
                'n_values': 3,
                'n_uniq_values': 3,
                'values_repartition': {
                    's1': 1,
                    's2': 1,
                    's3': 1
                }
            },
            {
                'col_name': 'col_2',
                'col_type': 'float64',
                'python_col_type': 'float',
                'n_values': 3,
                'n_uniq_values': 3,
                'mean': 15.9
            },
            {
                'col_name': 'col_3',
                'col_type': 'object',
                'python_col_type': 'str',
                'n_values': 3,
                'n_uniq_values': 2,
                'values_repartition': {
                    'M': 2,
                    'F': 1
                }
            }
        ]
        self.assertListEqual(self.parser.get_stats(), expected_list)
