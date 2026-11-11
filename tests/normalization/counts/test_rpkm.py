
from unittest import TestCase

import pandas as pd

from moonstone.normalization.counts.rpkm import (
    RPKMNormalization
)


class TestRPKMNormalization(TestCase):

    def setUp(self):
        self.raw_data = [
            [100, 200],
            [300, 0],
        ]
        self.column_names = ['Sample_1', 'Sample_2']
        self.index = ['Gene_1', 'Gene_2']

        self.raw_df = pd.DataFrame(
            self.raw_data, columns=self.column_names, index=self.index
        )

        self.feature_lengths = pd.Series(
            [1000, 2000], index=self.index
        )

    def test_normalize(self):
        # Verify RPKM calculations against manually computed expected values
        # Test automatic sequencing depth normalization
        # Check that zero counts remain zero
        # Verify preservation of sample and feature identifiers

        expected_data = [
            [250000.0, 1000000.0],
            [375000.0, 0.0],
        ]
        expected_df = pd.DataFrame(
            expected_data, columns=self.column_names, index=self.index
        )

        tested_normalization = RPKMNormalization(
            self.raw_df,
            feature_lengths=self.feature_lengths,
        )

        pd.testing.assert_frame_equal(
            tested_normalization.normalized_df, expected_df
        )

    def test_normalize_custom_depth(self):
        # Test custom sequencing depth normalization

        sequencing_depth = pd.Series(
            [1000, 2000], index=self.column_names
        )

        expected_data = [
            [100000.0, 100000.0],
            [150000.0, 0.0],
        ]
        expected_df = pd.DataFrame(
            expected_data, columns=self.column_names, index=self.index
        )

        tested_normalization = RPKMNormalization(
            self.raw_df,
            feature_lengths=self.feature_lengths,
            sequencing_depth=sequencing_depth,
        )

        pd.testing.assert_frame_equal(
            tested_normalization.normalized_df, expected_df
        )

    def test_invalid_feature_lengths(self):
        # Test missing, zero, or negative feature lengths

        invalid_lengths = [None, 0, -100, "invalid"]

        for value in invalid_lengths:
            with self.subTest(value=value):
                feature_lengths = self.feature_lengths.astype(object).copy()
                feature_lengths.loc["Gene_1"] = value

                tested_normalization = RPKMNormalization(
                    self.raw_df,
                    feature_lengths=feature_lengths,
                )

                with self.assertRaises(ValueError):
                    tested_normalization.normalized_df

    def test_missing_feature_length(self):
        # Test a feature missing from the feature lengths index

        feature_lengths = self.feature_lengths.drop("Gene_1")

        tested_normalization = RPKMNormalization(
            self.raw_df,
            feature_lengths=feature_lengths,
        )

        with self.assertRaises(ValueError):
            tested_normalization.normalized_df