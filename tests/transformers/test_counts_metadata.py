from unittest import TestCase

import pytest

import pandas as pd

from moonstone.transformers.mergers import (
    MergeCountsAndMetadata
)


class TestMergeCountsAndMetadata(TestCase):

    def test_check_column_names_True(self):
        tested_object_reads = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object_reads.columns.name = 'sample'
        tested_object_metadata = pd.DataFrame.from_dict(
            {
                'Sex': ['M', 'F', 'M', "F"],
                'AGE': [25, 65, 49, 50]
            },
            orient='index', columns=[1, 2, 3, 4])
        tested_object_metadata.columns.name = 'sample'
        tested_object = MergeCountsAndMetadata(tested_object_metadata, tested_object_reads)
        self.assertTrue(tested_object.check_column_names())

    def test_full_dataframe_with_different_sample_number(self):
        tested_object_reads = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object_reads.columns.name = 'sample'
        tested_object_metadata = pd.DataFrame.from_dict(
            {
                'Sex': ['M', 'F', 'M'],
                'AGE': [25, 65, 49]
            },
            orient='index', columns=[1, 2, 3])
        tested_object_metadata.columns.name = 'sample'
        tested_object = MergeCountsAndMetadata(tested_object_metadata, tested_object_reads)
        with pytest.raises(Exception):
            assert tested_object.full_dataframe

    def test_full_dataframe(self):
        tested_object_reads = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object_reads.columns.name = 'sample'
        tested_object_metadata = pd.DataFrame.from_dict(
            {
                'Sex': ['M', 'F', 'M', "F"],
                'AGE': [25, 65, 49, 50]
            },
            orient='index', columns=[1, 2, 3, 4])
        tested_object_metadata.columns.name = 'sample'
        tested_object = MergeCountsAndMetadata(tested_object_metadata, tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
                'Sex': ['M', 'F', 'M', "F"],
                'AGE': [25, 65, 49, 50],
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        expected_object.columns.name = 'sample'
        pd.testing.assert_frame_equal(tested_object.full_dataframe, expected_object)

    def test_full_dataframe_with_features_in_columns(self):
        tested_object_reads = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9]
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object_reads.columns.name = 'sample'
        tested_object_metadata = pd.DataFrame.from_dict(
            {
                'Sex': ['M', 'F', 'M', "F"],
                'AGE': [25, 65, 49, 50]
            },
            orient='index', columns=[1, 2, 3, 4])
        tested_object_metadata.columns.name = 'sample'
        tested_object = MergeCountsAndMetadata(tested_object_metadata, tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
                '1': ['M', 25, 3, 25],
                '2': ['F', 65, 2, 6],
                '3': ['M', 49, 1, 3],
                '4': ['F', 50, 0, 9]
            },
            orient='index', columns=['Sex', 'AGE', 'specie_1', 'specie_2'], dtype=str)
        expected_object.index.name = 'sample'
        pd.testing.assert_frame_equal(tested_object.full_dataframe_with_features_in_columns.astype(str),
                                      expected_object)
