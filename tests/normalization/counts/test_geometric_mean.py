from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.normalization.counts.geometric_mean import (
    GeometricMeanNormalization
)


class TestGeometricMeanNormalization(TestCase):

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
        tested_object = GeometricMeanNormalization(self.dummy_df)
        pd.testing.assert_frame_equal(tested_object.df, self.dummy_df)

    def test_non_zero_df(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=80)
        data = [
            [255, 26, 48, 75],
            [955, 222, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind).astype('float')
        pd.testing.assert_frame_equal(tested_object.non_zero_df(self.dummy_df), expected_result)

    def test_non_zero_df_threshold70(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=70, normalization_level=0)
        data = [
            [255, 26, 48, 75],
            [366, 46, 78, np.nan],
            [955, 222, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", 'Gen_2', "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind).astype('float')
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
        tested_object = GeometricMeanNormalization(self.dummy_df)

        data = [
            [5.541264, 3.258097, 3.871201, 4.317488],
            [6.861711, 5.402677, 3.828641, 4.174387],
            [4.488636, 3.988984, 4.976734, 3.367296]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)

        pd.testing.assert_frame_equal(tested_object.log_df(input_df), expected_result)

    def test_log_base_n_df(self):
        input_data = [
            [255, 26, 48, 75],
            [955, 222, 46, 65],
            [89, 54, 145, 29]
        ]
        input_column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        input_ind = ["Gen_1", "Gen_3", 'Gen_4']
        input_df = pd.DataFrame(input_data, columns=input_column_names, index=input_ind)
        tested_object = GeometricMeanNormalization(self.dummy_df, log_number=10)

        data = [
            [2.406540, 1.414973, 1.681241, 1.875061],
            [2.980003, 2.346353, 1.662758, 1.812913],
            [1.949390, 1.732394, 2.161368, 1.462398]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)

        pd.testing.assert_frame_equal(tested_object.log_df(input_df),
                                      expected_result)

    def test_removed_zero_df_None(self):
        data = [
            [255, 26, 48, 75],
            [366, 0, 78, 0],
            [955, 0, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", 'Gen_2', "Gen_3", 'Gen_4']
        dummy_df = pd.DataFrame(data, columns=column_names, index=ind)
        tested_object = GeometricMeanNormalization(dummy_df)
        expected_result = None
        self.assertEqual(tested_object.removed_zero_df, expected_result)

    def test_removed_zero_df(self):
        data = [
            [255, 26, 48, 75],
            [366, 0, 78, 0],
            [955, 0, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", 'Gen_2', "Gen_3", 'Gen_4']
        dummy_df = pd.DataFrame(data, columns=column_names, index=ind)
        tested_object = GeometricMeanNormalization(dummy_df)
        data = [
            [366, 0, 78, 0],
            [955, 0, 46, 65],
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ['Gen_2', "Gen_3"]
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)
        scaling_factors = tested_object.scaling_factors  # noqa
        pd.testing.assert_frame_equal(tested_object.removed_zero_df,
                                      expected_result)

    def test_calculating_and_substracting_mean_row(self):
        input_data = [
            [5.541264, 3.258097, 3.871201, 4.317488],
            [6.861711, 5.402677, 3.828641, 4.174387],
            [4.488636, 3.988984, 4.976734, 3.367296]
        ]
        input_column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        input_ind = ["Gen_1", "Gen_3", 'Gen_4']
        input_df = pd.DataFrame(input_data, columns=input_column_names, index=input_ind)
        tested_object = GeometricMeanNormalization(self.dummy_df)
        data = [
            [1.294251, -0.988916, -0.375811, 0.070475],
            [1.794857, 0.335823, -1.238213, -0.892467],
            [0.283224, -0.216428, 0.771321, -0.838117]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)

        pd.testing.assert_frame_equal(tested_object.calculating_and_substracting_mean_row(input_df).round(6),
                                      expected_result)

    def test_scaling_factor_zero_thresh_100(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=100)
        data = [3.648263, 0.805390, 0.686732, 0.432524]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.scaling_factors, expected_result)

    def test_scaling_factor_zero_thresh_80(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=80)
        data = [3.648263, 0.805390, 0.686732, 0.432524]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.scaling_factors, expected_result)

    def test_scaling_factor_zero_thresh_70(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=70)
        data = [3.495248, 0.612726, 0.699505, 0.432524]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.scaling_factors, expected_result)

    def test_scaling_factor_zero_thresh_70_more_zeros(self):
        data = [
            [255, 26, 48, 75],
            [366, 0, 78, 0],
            [955, 0, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", 'Gen_2', "Gen_3", 'Gen_4']
        more_zero_example = pd.DataFrame(data, columns=column_names, index=ind)
        tested_object = GeometricMeanNormalization(more_zero_example, zero_threshold=70)
        data = [3.648263, 0.588685, 0.686732, 0.458165]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.scaling_factors, expected_result)

    def test_scaling_factor_zero_thresh_80_more_zeros(self):
        data = [
            [255, 26, 48, 75],
            [366, 0, 78, 0],
            [955, 0, 46, 65],
            [89, 54, 145, 29]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", 'Gen_2', "Gen_3", 'Gen_4']
        more_zero_example = pd.DataFrame(data, columns=column_names, index=ind)
        tested_object = GeometricMeanNormalization(more_zero_example, zero_threshold=80)
        data = [2.487833, 0.588685, 1.424677, 0.752771]
        ind = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        expected_result = pd.Series(data, index=ind)
        pd.testing.assert_series_equal(tested_object.scaling_factors, expected_result)

    def test_normalized_df(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=100)
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

    def test_normalized_df_replace_0_to_1(self):
        tested_object = GeometricMeanNormalization(self.dummy_df, zero_threshold=100, replace_0_to_1=True)
        data = [
            [52.757471,	24.026807, 33.691852, 178.111683],
            [75.722488, 42.508965, 54.749259, 0.000000],
            [197.581903, 205.151963, 32.288025, 154.363458],
            [18.413392, 49.901829, 101.777469, 68.869851]
        ]
        column_names = ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4']
        ind = ["Gen_1", "Gen_2", "Gen_3", 'Gen_4']
        expected_result = pd.DataFrame(data, columns=column_names, index=ind)
        pd.testing.assert_frame_equal(tested_object.normalized_df, expected_result)

    def test_multiindex_normalization_first_element(self):
        tested_object = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                ('genus_3', 'specie_5'): {'1': 5, '2': 0, '3': 3, '4': 2},
                ('genus_3', 'specie_6'): {'1': 19, '2': 8, '3': 0, '4': 0},
            },
            orient='index')
        tested_object.columns.name = 'sample'
        tested_object.index.set_names(['genus', 'specie'], inplace=True)
        normalized_tested_object = GeometricMeanNormalization(tested_object, normalization_level='genus')
        expected_result = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 1.520270, '2': 1.454854, '3': 1.914414, '4': 0.000000},
                ('genus_1', 'specie_3'): {'1': 2.533784, '2': 5.819417, '3': 7.657658, '4': 10.050676},
                ('genus_2', 'specie_2'): {'1': 12.668919, '2': 4.364563, '3': 5.743243, '4': 12.922298},
                ('genus_2', 'specie_4'): {'1': 4.560811, '2': 13.093689, '3': 11.486487, '4': 4.307433},
                ('genus_3', 'specie_5'): {'1': 2.533784, '2': 0.000000, '3': 5.7432433, '4': 2.871622},
                ('genus_3', 'specie_6'): {'1': 9.628379, '2': 5.819417, '3': 0.000000, '4': 0.000000},
            },
            orient='index')
        expected_result.columns.name = 'sample'
        expected_result.index.set_names(['genus', 'specie'], inplace=True)
        pd.testing.assert_frame_equal(normalized_tested_object.normalized_df, expected_result)

    def test_multiindex_normalization_second_element(self):
        tested_object = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                ('genus_3', 'specie_5'): {'1': 5, '2': 0, '3': 3, '4': 2},
                ('genus_3', 'specie_6'): {'1': 19, '2': 8, '3': 0, '4': 0},
            },
            orient='index')
        tested_object.columns.name = 'sample'
        tested_object.index.set_names(['genus', 'specie'], inplace=True)
        normalized_tested_object = GeometricMeanNormalization(tested_object, normalization_level=1)
        expected_result = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 2.449490, '2': 1.446254, '3': 1.446254, '4': 0.000000},
                ('genus_1', 'specie_3'): {'1': 4.082483, '2': 5.785015, '3': 5.785015, '4': 6.204679},
                ('genus_2', 'specie_2'): {'1': 20.412415, '2': 4.338761, '3': 4.338761, '4': 7.977444},
                ('genus_2', 'specie_4'): {'1': 7.348469, '2': 13.016284, '3': 8.677523, '4': 2.659148},
                ('genus_3', 'specie_5'): {'1': 4.082483, '2': 0.000000, '3': 4.338761, '4': 1.772765},
                ('genus_3', 'specie_6'): {'1': 15.513435, '2': 5.785015, '3': 0.000000, '4': 0.000000},
            },
            orient='index')
        expected_result.columns.name = 'sample'
        expected_result.index.set_names(['genus', 'specie'], inplace=True)
        pd.testing.assert_frame_equal(normalized_tested_object.normalized_df, expected_result)
