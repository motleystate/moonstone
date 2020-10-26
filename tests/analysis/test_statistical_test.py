from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.analysis.statistical_test import (
    mann_whitney_u_group
)


class TestShannonIndex(TestCase):

    def test_mann_whitney_u_group(self):

        test_df = pd.Series({
            'sample1': 0.72,
            'sample2': 0.98,
            'sample3': 0.00,
            'sample4': 0.00,
            'sample5': 1.99
        })
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )

        expected_df = pd.DataFrame(
            [
                [np.nan, 0.374259],
                [0.374259, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )

        matrix = mann_whitney_u_group(test_df, metadata_df['group'])
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_group_wrong_metadata(self):
        test_df = pd.Series({
            'sample1': 0.72,
            'sample2': 0.98,
            'sample3': 0.00,
            'sample4': 0.00,
            'sample5': 1.99
        })
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1]
            ],
            columns=['sex', 'group'],
            index=['sample6', 'sample7'],
        )

        expected_df = pd.DataFrame(
            [
                [np.nan, 0],
                [0, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )

        matrix = mann_whitney_u_group(test_df, metadata_df['group'])
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)
