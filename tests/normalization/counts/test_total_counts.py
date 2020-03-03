from unittest import TestCase

import pandas as pd

from moonstone.normalization.counts.total_counts import (
    TotalCountsNormalization
)


class TestTotalCountsNormalization(TestCase):

    def setUp(self):
        self.raw_data = [
            [199, 1, 48, 75],
            [0, 24, 1, 0],
            [1, 25, 1, 25],
        ]
        self.column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        self.index = ["Gen_1", 'Gen_2', "Gen_3"]
        self.raw_df = pd.DataFrame(self.raw_data, columns=self.column_names, index=self.index)

    def test_normalize(self):
        tested_data = [
            [99.5, 2.0, 96.0, 75.0],
            [0.0, 48.0, 2.0, 0.0],
            [1.0, 50.0, 2.0, 25.0],
        ]
        expected_df = pd.DataFrame(tested_data, columns=self.column_names, index=self.index)
        tested_normalization = TotalCountsNormalization(self.raw_df)
        pd.testing.assert_frame_equal(tested_normalization.normalized_df, expected_df)
        # Check scaling factors
        expected_scaling_factors = pd.Series(
            [2, 0.5, 0.5, 1], index=self.column_names
        )
        pd.testing.assert_series_equal(tested_normalization.scaling_factors, expected_scaling_factors)
