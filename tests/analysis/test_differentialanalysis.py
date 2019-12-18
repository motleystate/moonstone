from unittest import TestCase
import pytest

import pandas as pd
import numpy as np

from moonstone.analysis.differentialanalysis import (
    DifferentialAnalysis
)


class TestDifferentialAnalysis(TestCase):
    def setUp(self):
        self.tested_object_reads = pd.DataFrame.from_dict(
            {
                'specie_1': [300.0, 245.0, 11.0, 154.0, 239.0, 789.0, 81.0, 58.0, 76.0, 28.0],
                'specie_2': [3886.0, 1184.0, 1250.0, 6347.0, 1208.0, 3403.0, 1522.0, 1378.0, 2333.0, 2749.0],
                'specie_3': [0.0, 0.0, 0.0, 0.0, 1961.0, 0.0, 0.0, 30.0, 0.0, 0.0],
                'specie_4': [2362.0, 1858.0, 23.0, 1300.0, 515.0, 17.0, 39.0, 1111.0, 549.0, 54.0],
                'specie_5': [0.0, 295.0, 472.0, 0.0, 23.0, 740.0, 136.0, 39.0, 1.0, 0.0]
            },
            orient='index', columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        self.tested_object_reads.columns.name = 'sample'
        self.tested_object_metadata = pd.DataFrame.from_dict(
            {
                'SEX': ['M', 'M', 'F', "M", 'M', 'M', 'F', 'F', 'F', 'F'],
                'Season': ['Fall', 'Fall', 'Spring', 'Summer', 'Summer', 'Spring', 'Spring', 'Fall', 'Summer', 'Spring']
            },
            orient='index', columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        self.tested_object_metadata.columns.name = 'sample'
        self.dicotomic_feature = 'SEX'
        self.dicotomic_features = ['SEX']
        self.multiple_option_feature = 'Season'
        self.multiple_option_features = ['Season']
        self.significance_level = 0.05

    def test_number_columns_to_skip(self):
        tested_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = 2
        self.assertEqual(tested_object.number_columns_to_skip, expected_object)

    def test_test_dichotomic_features_t_test(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
                0: ['SEX', 'specie_1', 2.580936, 0.032569, 64219.3, 925.7],
                1: ['SEX', 'specie_2', 1.353146, 0.213002, 4611894.3, 432948.3],
                2: ['SEX', 'specie_3', 0.984586, 0.353664, 769104.2, 180.0],
                3: ['SEX', 'specie_4', 1.788817, 0.111441, 915345.3, 227463.2],
                4: ['SEX', 'specie_5', 0.485682, 0.640215, 102820.3, 39705.3],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value', 'variance_group1',
                                     'variance_group2'])
        pd.testing.assert_frame_equal(test_object.test_dichotomic_features(self.dicotomic_feature, 't_test'),
                                      expected_object)

    def test_test_dichotomic_features_wilcoxon_rank_test(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['SEX', 'specie_1', 2.611165, 0.009023, 64219.3, 925.7],
               1: ['SEX', 'specie_2', 0.522233, 0.601508, 4611894.3, 432948.3],
               2: ['SEX', 'specie_3', 0.104447, 0.916815, 769104.2, 180.0],
               3: ['SEX', 'specie_4', 1.148913, 0.250592, 915345.3,	227463.2],
               4: ['SEX', 'specie_5', -0.104447, 0.916815, 102820.3, 39705.3],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value', 'variance_group1',
                                     'variance_group2'])
        pd.testing.assert_frame_equal(test_object.test_dichotomic_features(self.dicotomic_feature,
                                                                           'wilcoxon_rank_test'), expected_object)

    def test_test_multiple_features_one_way_anova(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['Season', 'specie_1', 0.064622, 0.937974],
               1: ['Season', 'specie_2', 0.401588, 0.683744],
               2: ['Season', 'specie_3', 1.208941, 0.354003],
               3: ['Season', 'specie_4', 15.370573, 0.002748],
               4: ['Season', 'specie_5', 1.817575, 0.231335],

            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value'])
        pd.testing.assert_frame_equal(test_object.test_multiple_features(self.multiple_option_feature, 'one_way_anova'),
                                      expected_object)

    def test_test_multiple_features_kruskal_test(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['Season', 'specie_1', 0.890909, 0.640533],
               1: ['Season', 'specie_2', 0.336364, 0.845200],
               2: ['Season', 'specie_3', 1.518519, 0.468013],
               3: ['Season', 'specie_4', 7.436364, 0.024278],
               4: ['Season', 'specie_5', 2.142857, 0.342519],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value'])
        pd.testing.assert_frame_equal(test_object.test_multiple_features(self.multiple_option_feature, 'kruskal_test'),
                                      expected_object)

    def test_test_default(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        with pytest.raises(Exception):
            assert test_object.test_default()

    def test_corrected_p_values(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        tested_object = pd.DataFrame.from_dict(
            {
               0: ['Season', 'specie_1', 0.890909, 0.640533],
               1: ['Season', 'specie_2', 0.336364, 0.845200],
               2: ['Season', 'specie_3', 1.518519, 0.468013],
               3: ['Season', 'specie_4', 7.436364, 0.024278],
               4: ['Season', 'specie_5', 2.142857, 0.342519],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value'])
        expected_object = np.array([0.8006663, 0.8452, 0.7800217, 0.12139, 0.7800217])
        np.testing.assert_almost_equal(test_object.corrected_p_values(tested_object['p-value'], 'fdr_bh'),
                                       expected_object)

    def test_differential_analysis_by_feature(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['Season', 'specie_1', 0.890909, 0.640533, 0.800666],
               1: ['Season', 'specie_2', 0.336364, 0.845200, 0.845200],
               2: ['Season', 'specie_3', 1.518519, 0.468013, 0.780022],
               3: ['Season', 'specie_4', 7.436364, 0.024278, 0.121390],
               4: ['Season', 'specie_5', 2.142857, 0.342519, 0.780022],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value', 'corrected_p-value'])
        pd.testing.assert_frame_equal(test_object.differential_analysis_by_feature(self.multiple_option_features,
                                                                                   'multiple_features',
                                                                                   'kruskal_test',
                                                                                   'fdr_bh'),
                                      expected_object)
