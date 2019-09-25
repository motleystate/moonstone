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
                'col_3': ['M', 'F', 'M'],
                'col_4': ['T', 'F', 'F']
            }
        )
        assert_frame_equal(self.parser.dataframe, expected_df)


class TestMetadataParserBuildStats(TestCase):

    def setUp(self):
        #  Instantiate a parser with no file to call protected methods
        self.parser = MetadataParser("no_file")
        self.parser._dataframe = pd.DataFrame(
            {
                'test_int': [12, 41, 15],
                'test_float': [13.3, 15.3, 19.1],
                'test_str': ['M', 'F', 'M']
            }
        )

    def test_build_stats_str(self):
        selected_col = 'test_str'
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
        print(self.parser._build_stats(selected_col))
        self.assertDictEqual(self.parser._build_stats(selected_col), expected_dict)

    def test_build_stats_int(self):
        selected_col = 'test_int'
        expected_dict = {
            'col_name': 'test_int',
            'col_type': 'int64',
            'python_col_type': 'int',
            'n_values': 3,
            'n_uniq_values': 3,
            'mean': 22.67
        }
        self.assertDictEqual(self.parser._build_stats(selected_col), expected_dict)

    def test_build_stats_float(self):
        selected_col = 'test_float'
        expected_dict = {
            'col_name': 'test_float',
            'col_type': 'float64',
            'python_col_type': 'float',
            'n_values': 3,
            'n_uniq_values': 3,
            'mean': 15.9
        }
        self.assertDictEqual(self.parser._build_stats(selected_col), expected_dict)
