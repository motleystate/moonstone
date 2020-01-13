import os
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from moonstone.parsers.metadata import MetadataParser


class TestMetadataParser(TestCase):

    def setUp(self):
        self.metadata_file = os.path.join(os.path.dirname(__file__), "data/metadata/metadata.tsv")
        self.metadata_file_no_header = os.path.join(os.path.dirname(__file__), "data/metadata/metadata_noheader.tsv")
        self.metadata_file_dirty = os.path.join(os.path.dirname(__file__), "data/metadata/dirty_metadata.tsv")

    def test_parse_file(self):
        expected_df = pd.DataFrame(
            {
                'col_1': ['s1', 's2', 's3'],
                'col_2': [13.3, 15.3, 19.1],
                'col_3': ['M', 'F', 'M']
            }
        )
        parser = MetadataParser(self.metadata_file)
        assert_frame_equal(parser.dataframe, expected_df)

    def test_get_stats_headers(self):
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
                'mean': 15.9,
                'values_repartition': {
                    13.3: 1,
                    15.3: 1,
                    19.1: 1
                }
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
        parser = MetadataParser(self.metadata_file)
        self.assertListEqual(parser.get_stats(), expected_list)

    def test_get_stats_no_header(self):
        expected_list = [
            {
                'col_name': 0,
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
                'col_name': 1,
                'col_type': 'float64',
                'python_col_type': 'float',
                'n_values': 3,
                'n_uniq_values': 3,
                'mean': 15.9,
                'values_repartition': {
                    13.3: 1,
                    15.3: 1,
                    19.1: 1
                }
            },
            {
                'col_name': 2,
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
        parser = MetadataParser(self.metadata_file_no_header, no_header=True)
        self.assertListEqual(parser.get_stats(), expected_list)

    def test_parse_file_force_dtype(self):
        expected_df = pd.DataFrame(
            {
                'col_1': ['s1', 's2', 's3'],
                'col_2': ['13.3', '15.3', '19.1'],
                'col_3': ['M', 'F', 'M']
            }
        )
        parsing_options = {
            'dtype': {
                'col_2': 'object'
            }
        }
        parser = MetadataParser(self.metadata_file, parsing_options=parsing_options)
        assert_frame_equal(parser.dataframe, expected_df)

    def test_parse_dirty_metadata_and_clean(self):
        expected_df = pd.DataFrame(
            {
                'sample': ['s1', 's2', 's3', 's4'],
                'age': [29, 48, 36, 25],
            }
        )
        cleaning_operations = {
            'samples': [
                ('to_slug', {}),
                ('rename', {'new_name': 'sample'})
            ]
        }
        parser = MetadataParser(self.metadata_file_dirty, sep=",", cleaning_operations=cleaning_operations)
        assert_frame_equal(parser.dataframe, expected_df)
