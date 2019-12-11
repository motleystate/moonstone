from unittest import TestCase

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
        self.dicotomic_features = ['SEX']
        self.multiple_option_features = ['Season']
        self.significance_level = 0.05

    def test_full_table(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               '1': ['M', 'Fall', 300.0, 3886.0, 0.0, 2362.0, 0.0],
                '2': ['M', 'Fall', 245.0, 1184.0, 0.0, 1858.0, 295.0],
                '3': ['F', 'Spring', 11.0, 1250.0, 0.0, 23.0, 472.0],
                '4': ['M', 'Summer', 154.0, 6347.0, 0.0, 1300.0, 0.0],
                '5': ['M', 'Summer', 239.0, 1208.0, 1961.0, 515.0, 23.0],
                '6': ['M', 'Spring', 789.0, 3403.0, 0.0, 17.0, 740.0],
                '7': ['F', 'Spring', 81.0, 1522.0, 0.0, 39.0, 136.0],
                '8': ['F', 'Fall', 58.0, 1378.0, 30.0, 1111.0, 39.0],
                '9': ['F', 'Summer', 76.0, 2333.0, 0.0, 549.0, 1.0],
                '10': ['F', 'Spring', 28.0, 2749.0, 0.0, 54.0, 0.0],
            },
            orient='index', columns=['SEX', 'Season', 'specie_1', 'specie_2', 'specie_3', 'specie_4', 'specie_5'], dtype=str)
        expected_object.index.name = 'sample'
        pd.testing.assert_frame_equal(test_object.full_table.astype(str), expected_object)

    def test_number_columns_to_skip(self):
        tested_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = 2
        self.assertEqual(tested_object.number_columns_to_skip, expected_object)

    def test_t_test(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
                0: ['SEX', 'specie_1', 2.580936, 0.032569, 64219.3, 925.7],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value', 'variance_group1', 'variance_group2'])
        pd.testing.assert_frame_equal(test_object.t_test(self.dicotomic_features, self.significance_level), expected_object)

    def test_wilcoxon_rank_test(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['SEX', 'specie_1', 2.611165, 0.009023, 64219.3, 925.7],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value', 'variance_group1', 'variance_group2'])
        pd.testing.assert_frame_equal(test_object.wilcoxon_rank_test(self.dicotomic_features, self.significance_level), 
                                      expected_object)

    def test_one_way_anova(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['Season', 'specie_4', 15.370573, 0.002748],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value'])
        pd.testing.assert_frame_equal(test_object.one_way_anova(self.multiple_option_features, self.significance_level),
                                      expected_object)

    def test_kruskal_test(self):
        test_object = DifferentialAnalysis(self.tested_object_metadata, self.tested_object_reads)
        expected_object = pd.DataFrame.from_dict(
            {
               0: ['Season', 'specie_4', 7.436364, 0.024278],
            },
            orient='index', columns=['features', 'taxons', 'static_value', 'p-value'])
        pd.testing.assert_frame_equal(test_object.kruskal_test(self.multiple_option_features, self.significance_level),
                                      expected_object)
