from unittest import TestCase

import pandas as pd

from moonstone.normalization.DESeq2_normalization import (
    DESeq2_normalization
)


class TestSequence(TestCase):

    def setUp(self):
        data = [
            [255, 26, 48, 75],
            [366, 46, 78, 0],
            [955, 222, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", 'Gen_2', "Gen_3", 'Gen_4']
        self.dummy_df = pd.DataFrame(data, columns=column_names, index=ind)

    def test_check_format(self):
        tested_object = DESeq2_normalization(self.dummy_df)
        pd.testing.assert_frame_equal(tested_object.df, self.dummy_df)

    def test_non_zero_df(self):
        tested_object = DESeq2_normalization(self.dummy_df)
        data = [
            [255, 26, 48, 75],
            [955, 222, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)
        pd.testing.assert_frame_equal(tested_object.non_zero_df(self.dummy_df), expected_result)

    def test_log_df(self):
        input_data = [
            [255, 26, 48, 75],
            [955, 222, 46, 65],
            [89, 54, 145, 29]
        ]
        input_column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        input_ind = ["Gen_1", "Gen_3", 'Gen_4']
        input_df = pd.DataFrame(input_data, columns=input_column_names, index=input_ind)
        tested_object = DESeq2_normalization(self.dummy_df)

        data = [
            [5.541264, 3.258097, 3.871201, 4.317488],
            [6.861711, 5.402677, 3.828641, 4.174387],
            [4.488636, 3.988984, 4.976734, 3.367296]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)

        pd.testing.assert_frame_equal(tested_object.log_df(input_df), expected_result)

    def test_calculating_and_substracting_mean_row(self):
        input_data = [
            [5.541264, 3.258097, 3.871201, 4.317488],
            [6.861711, 5.402677, 3.828641, 4.174387],
            [4.488636, 3.988984, 4.976734, 3.367296]
        ]
        input_column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        input_ind = ["Gen_1", "Gen_3", 'Gen_4']
        input_df = pd.DataFrame(input_data, columns=input_column_names, index=input_ind)
        tested_object = DESeq2_normalization(self.dummy_df)
        data = [
            [1.294251, -0.988916, -0.375811, 0.070475],
            [1.794857, 0.335823, -1.238213, -0.892467],
            [0.283224, -0.216428, 0.771321, -0.838117]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)

        pd.testing.assert_frame_equal(tested_object.calculating_and_substracting_mean_row(input_df).round(6), expected_result)

    #not changed yet
    def test_sample_log_median(self):
        tested_object = DESeq2_normalization(self.dummy_df)
        data = [1.294251, -0.216428, -0.375811, -0.838117]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.sample_log_median(), expected_result)

    def test_scaling_factor(self):
        tested_object = DESeq2_normalization(self.dummy_df)
        data = [3.648263, 0.805390, 0.686732, 0.432524]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.scaling_factor, expected_result)

    def test_normalized_df(self):
        tested_object = DESeq2_normalization(self.dummy_df)
        data = [
            [69.896271,	32.282490, 69.896271, 173.400644],
            [100.321707, 57.115175, 113.581441, 0.000000],
            [261.768388, 275.642802, 66.983926, 150.280558],
            [24.395169, 67.048249, 211.144986, 67.048249]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_2", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)
        pd.testing.assert_frame_equal(tested_object.normalized_df, expected_result)
