from unittest import TestCase

import pandas as pd

from moonstone.normalization.counts.random_selection import (
    RandomSelection
)


class TestRandomSelection(TestCase):

    def setUp(self):
        self.raw_data = [
            [199, 1, 48, 75],
            [0, 24, 1, 0],
            [1, 25, 1, 25],
        ]
        self.column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        self.index = ["Gen_1", 'Gen_2', "Gen_3"]
        self.raw_df = pd.DataFrame(self.raw_data, columns=self.column_names, index=self.index)

    def test_normalize_default_threshold(self):
        expected_data = [
            [50, 1, 48, 40],
            [0, 24, 1, 0],
            [0, 25, 1, 10],
        ]
        expected_df = pd.DataFrame(expected_data, columns=self.column_names, index=self.index)
        tested_normalization = RandomSelection(self.raw_df, random_seed=2935)
        pd.testing.assert_frame_equal(tested_normalization.normalized_df, expected_df)

    def test_normalize_threshold_20(self):
        expected_data = [
            [20, 0, 20, 16],
            [0, 9, 0, 0],
            [0, 11, 0, 4],
        ]
        expected_df = pd.DataFrame(expected_data, columns=self.column_names, index=self.index)
        tested_normalization = RandomSelection(self.raw_df, threshold=20, random_seed=2935)
        pd.testing.assert_frame_equal(tested_normalization.normalized_df, expected_df)

    def test_normalize_threshold_100(self):
        # Define what we expect here, do we want to filter out instead below threshold?
        expected_data = [
            [100, 1, 48, 75],
            [0, 24, 1, 0],
            [0, 25, 1, 25],
        ]
        expected_df = pd.DataFrame(expected_data, columns=self.column_names, index=self.index)
        tested_normalization = RandomSelection(self.raw_df, threshold=100, random_seed=2935)
        pd.testing.assert_frame_equal(tested_normalization.normalized_df, expected_df)
