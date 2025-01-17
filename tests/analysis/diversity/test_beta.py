from unittest import TestCase

from io import StringIO
import numpy as np
import pandas as pd
from skbio import TreeNode

from moonstone.analysis.diversity.beta import (
    BrayCurtis, WeightedUniFrac, UnweightedUniFrac
)


class TestBrayCurtis(TestCase):

    def setUp(self):
        self.tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0],
                'species2': [1, 0, 2],
                'species3': [0, 0, 0],
                'species4': [0, 3, 0]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )

    def test_compute_beta_diversity_df(self):
        tested_object_instance = BrayCurtis(self.tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'sample1': [0, 0.333, 0.714],
                'sample2': [0.333, 0, 1.0],
                'sample3': [0.714, 1.0, 0],
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )
        pd.testing.assert_frame_equal(
            tested_object_instance.beta_diversity_df, expected_object,
            rtol=0.01
            # check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
        )

    def test_compute_beta_diversity_series(self):

        tested_object_instance = BrayCurtis(self.tested_object)
        multi_index = pd.MultiIndex.from_tuples(
            [
                ('sample1', 'sample2'),
                ('sample1', 'sample3'),
                ('sample2', 'sample3'),
            ]
        )
        expected_object = pd.Series(
            [0.333, 0.714, 1.0], index=multi_index, name=BrayCurtis._DIVERSITY_INDEXES_NAME
        )
        # Two ways of retrieving the series
        pd.testing.assert_series_equal(
            tested_object_instance.beta_diversity_series, expected_object,
            rtol=0.01
        )
        pd.testing.assert_series_equal(
            tested_object_instance.diversity_indexes, expected_object,
            rtol=0.01
        )

    def test_run_statistical_test_groups_with_NaN(self):
        tested_object_instance = BrayCurtis(self.tested_object)

        tested_df = pd.DataFrame.from_dict(
            {
                'samples1': [1.88, 'B'],
                'samples2': [2.43, 'C'],
                'samples3': [7.98, 'C'],
                'samples4': [2.89, 'A'],
                'samples5': [0.07, 'B'],
                'samples6': [4.77, 'B'],
                'samples7': [9.76, 'A'],
                'samples8': [9.40, 'A'],
                'samples9': [2.34, 'B'],
                'samples10': [5.67, 'A'],
                'samples11': [1.26, 'A'],
                'samples12': [2.31, 'C'],
                'samples13': [1.19, 'B'],
                'samples14': [9.35, 'A'],
                'samples15': [7.89, 'A'],
                'samples16': [4.65, 'C'],
                'samples17': [8.90, 'D'],  # only 1 sample from group D < 5 required to do ttest-independence
                'samples18': [2.33, 'C'],
                'samples19': [1.34, 'B'],
                'samples20': [6.87, 'C']
            },
            orient='index', columns=['beta_index', 'Group']
        )

        expected_object = pd.Series(
            {
                ('A', 'B'): 0.03220171408367315,
                ('A', 'C'): 0.21475494241030912,
                ('B', 'C'): 0.10740702992087024,
                ('A', 'D'): np.nan,
                ('B', 'D'): np.nan,
                ('C', 'D'): np.nan
            }
        )
        expected_object.index.names = ["Group1", "Group2"]

        pval = tested_object_instance._run_statistical_test_groups(
            tested_df, 'Group', stats_test='ttest_independence', correction_method='fdr_bh',
            structure_pval='series', sym=False
            )

        pd.testing.assert_series_equal(pval, expected_object)

        # and checking structure_pval='dataframe'
        expected_object = pd.DataFrame.from_dict({
            'B': [0.03220171408367315, np.nan, np.nan],
            'C': [0.21475494241030912, 0.10740702992087024, np.nan],
            'D': [np.nan, np.nan, np.nan]
        })
        expected_object.index = ["A", "B", "C"]

        pval = tested_object_instance._run_statistical_test_groups(
            tested_df, 'Group', stats_test='ttest_independence', correction_method='fdr_bh',
            structure_pval='dataframe', sym=False
            )

        pd.testing.assert_frame_equal(pval, expected_object)

    def test_get_grouped_df_series(self):
        metadata_ser = pd.Series(
            {
                'sample1': 'M',
                'sample2': 'F',
                'sample3': 'F',
            },
            name='sex'
        )
        tested_object_instance = BrayCurtis(self.tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'sample2-sample3': [1.0, 'F'],
            },
            orient='index', columns=['beta_index', 'sex']
        )
        output = tested_object_instance._get_grouped_df_series(metadata_ser)
        pd.testing.assert_frame_equal(
            output, expected_object,
        )

    def test_get_grouped_df_dataframe(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 5, 3, 1, 0],
                'species2': [1, 0, 2, 2, 1, 2, 2],
                'species3': [0, 0, 0, 0, 1, 0, 0],
                'species4': [0, 3, 0, 0, 2, 4, 1]
            },
            orient='index',
            columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7']
        )

        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M', 'A', 'M - A'],
                'sample2': ['F', 'B', 'F - B'],
                'sample3': ['F', 'B', 'F - B'],
                'sample4': ['F', 'A', 'F - A'],
                'sample5': ['M', 'B', 'M - B'],
                'sample6': ['M', 'A', 'M - A'],
                'sample7': ['F', 'A', 'F - A']
            },
            orient='index', columns=['sex', 'Group', 'sex_Group']
        )
        tested_object_instance = BrayCurtis(tested_object)

        expected_object = pd.DataFrame.from_dict(
            {
                'sample1-sample6': [2/3, 'M - A', 'M', 'A'],
                'sample2-sample3': [1.0, 'F - B', 'F', 'B'],
                'sample4-sample7': [0.6, 'F - A', 'F', 'A']
            },
            orient='index', columns=['beta_index', 'sex_Group', 'sex', 'Group']
        )

        output = tested_object_instance._get_grouped_df_dataframe(metadata_df)
        pd.testing.assert_frame_equal(
            output, expected_object,
        )

    def test_analyse_grouped_df(self):
        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M'],
                'sample2': ['F'],
                'sample3': ['F'],
            },
            orient='index', columns=['sex']
        )
        tested_object_instance = BrayCurtis(self.tested_object)

        expected_object = pd.DataFrame.from_dict(
            {
                'sample2-sample3': [1.0, 'F'],
            },
            orient='index', columns=['beta_index', 'sex']
        )
        output = tested_object_instance.analyse_groups(metadata_df, 'sex', show=False, show_pval=False)
        pd.testing.assert_frame_equal(
            output['data'], expected_object,
        )

    def test_analyse_grouped_df_with_group_col2(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 5, 3, 1, 0],
                'species2': [1, 0, 2, 2, 1, 2, 2],
                'species3': [0, 0, 0, 0, 1, 0, 0],
                'species4': [0, 3, 0, 0, 2, 4, 1]
            },
            orient='index',
            columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7']
        )

        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M', 'A'],
                'sample2': ['F', 'B'],
                'sample3': ['F', 'B'],
                'sample4': ['F', 'A'],
                'sample5': ['M', 'B'],
                'sample6': ['M', 'A'],
                'sample7': ['F', 'A']
            },
            orient='index', columns=['sex', 'Group']
        )
        tested_object_instance = BrayCurtis(tested_object)

        expected_object = pd.DataFrame.from_dict(
            {
                'sample1-sample6': [2/3, 'M - A', 'M', 'A'],
                'sample2-sample3': [1.0, 'F - B', 'F', 'B'],
                'sample4-sample7': [0.6, 'F - A', 'F', 'A']
            },
            orient='index', columns=['beta_index', 'sex_Group', 'sex', 'Group']
        )

        output = tested_object_instance.analyse_groups(
            metadata_df, group_col='sex', group_col2='Group', show=False, show_pval=False
            )
        pd.testing.assert_frame_equal(
            output["data"], expected_object,
        )

    def test_unit_scale_with_scale_between0and1(self):
        tested_object_instance = BrayCurtis(self.tested_object)

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [-30],
                'species2': [-1],
                'species3': [-0.2],
                'species4': [0],
                'species5': [0.5],
                'species6': [1],
                'species7': [198]
            },
            orient='index',
            columns=['PC1']
        )

        expected_object = pd.DataFrame.from_dict(
            {
                'species1': [-0.06],
                'species2': [-0.002],
                'species3': [-0.0004],
                'species4': [0],
                'species5': [0.001],
                'species6': [0.002],
                'species7': [0.396]
            },
            orient='index',
            columns=['PC1']
        )
        pd.testing.assert_frame_equal(
            tested_object_instance._unit_scale(tested_object, 228, 0.456),
            expected_object
        )

    def test_unit_scale_with_scale_above1(self):
        tested_object_instance = BrayCurtis(self.tested_object)

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [-30],
                'species2': [-1],
                'species3': [-0.2],
                'species4': [0],
                'species5': [0.5],
                'species6': [1],
                'species7': [198]
            },
            orient='index',
            columns=['PC1']
        )

        expected_object = pd.DataFrame.from_dict(
            {
                'species1': [-3],
                'species2': [-0.1],
                'species3': [-0.02],
                'species4': [0],
                'species5': [0.05],
                'species6': [0.1],
                'species7': [19.80]
            },
            orient='index',
            columns=['PC1']
        )
        pd.testing.assert_frame_equal(
            tested_object_instance._unit_scale(tested_object, 228, 22.8),
            expected_object
        )

    def test_scale_biplot(self):
        tested_object_instance = BrayCurtis(self.tested_object)

        tested_features_df = pd.DataFrame.from_dict(
            {
                ('genusA', 'species1'): [-30, 0],
                ('genusA', 'species2'): [-1, 89],
                ('genusB', 'species3'): [-0.2, -5.7],
                ('genusB', 'species4'): [0, 0],
                ('genusB', 'species5'): [0.5, 0.32],
                ('genusC', 'species6'): [1, -0.6],
                ('genusC', 'species7'): [198, 1]
            },
            orient='index',
            columns=['PC1', 'PC2']
        )
        tested_features_df.index = pd.MultiIndex.from_tuples(tested_features_df.index, names=('genus', 'species'))

        tested_samples_df = pd.DataFrame.from_dict(
            {
                'sample1': [-6, 3.4],
                'sample2': [2, 0.9],
                'sample3': [5, -2.4],
                'sample4': [13, -0.1],
                'sample5': [-0.5, 0],
                'sample6': [-62, -25],
                'sample7': [0, 0.3]
            },
            orient='index',
            columns=['PC1', 'PC2']
        )

        expected_object = pd.DataFrame.from_dict(
            {
                ('genusA', 'species1'): [-62, -23.290602],
                ('genusA', 'species2'): [-52.460526, 3.4],
                ('genusB', 'species3'): [-52.197368, -25],
                ('genusB', 'species4'): [-52.131579, -23.290602],
                ('genusB', 'species5'): [-51.967105, -23.194636],
                ('genusC', 'species6'): [-51.802632, -23.470539],
                ('genusC', 'species7'): [13, -22.990707]
            },
            orient='index',
            columns=['PC1', 'PC2']
        )
        expected_object.index = pd.MultiIndex.from_tuples(expected_object.index, names=('genus', 'species'))

        pd.testing.assert_frame_equal(
            tested_object_instance._scale_biplot(tested_features_df, tested_samples_df),
            expected_object
        )


class TestWeightedUniFrac(TestCase):

    def setUp(self):
        self.tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0],
                'species2': [1, 0, 2],
                'species3': [0, 0, 0],
                'species4': [0, 3, 0]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )

    def test_compute_beta_diversity(self):
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species4:0.5,species5:0.5):1.0)root;'))
        tested_object_instance = WeightedUniFrac(self.tested_object, tree)
        expected_object = pd.DataFrame.from_dict(
            {
                'sample1': [0, 1.285714, 0.400000],
                'sample2': [1.285714, 0, 1.571429],
                'sample3': [0.400000, 1.571429, 0],
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )
        pd.testing.assert_frame_equal(
            tested_object_instance.beta_diversity_df, expected_object,
        )

    def test_compute_beta_diversity_force_computation(self):
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species4:0.5,species5:0.5):1.0)root;'))
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 2, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4],
                'species6': [3, 5, 2, 2, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = WeightedUniFrac(tested_object, tree, force_computation=True)

        expected_object_instance = WeightedUniFrac(tested_object.iloc[:4], tree)
        expected_object = expected_object_instance.beta_diversity_df

        with self.assertLogs('moonstone.analysis.diversity.beta', level='WARNING') as log:
            tested_object = tested_object_instance.beta_diversity_df
            pd.testing.assert_frame_equal(tested_object, expected_object)

            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.analysis.diversity.beta:INCOMPLETE TREE: missing ['species6'].\n\
Computation of the Weighted UniFrac diversity using only the OTU IDs present in the Tree.", log.output)


class TestUnweightedUniFrac(TestCase):

    def setUp(self):
        self.tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0],
                'species2': [1, 0, 2],
                'species3': [0, 0, 0],
                'species4': [0, 3, 0]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )

    def test_compute_beta_diversity(self):
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species4:0.5,species5:0.5):1.0)root;'))
        tested_object_instance = UnweightedUniFrac(self.tested_object, tree)
        expected_object = pd.DataFrame.from_dict(
            {
                'sample1': [0, 0.538462, 0.142857],
                'sample2': [0.538462, 0, 0.615385],
                'sample3': [0.142857, 0.615385, 0],
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )
        pd.testing.assert_frame_equal(
            tested_object_instance.beta_diversity_df, expected_object,
        )

    def test_compute_beta_diversity_force_computation(self):
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species4:0.5,species5:0.5):1.0)root;'))
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 2, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4],
                'species6': [3, 5, 2, 2, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = UnweightedUniFrac(tested_object, tree, force_computation=True)

        expected_object_instance = UnweightedUniFrac(tested_object.iloc[:4], tree)
        expected_object = expected_object_instance.beta_diversity_df

        with self.assertLogs('moonstone.analysis.diversity.beta', level='WARNING') as log:
            tested_results = tested_object_instance.beta_diversity_df
            pd.testing.assert_frame_equal(tested_results, expected_object)

            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.analysis.diversity.beta:INCOMPLETE TREE: missing ['species6'].\n\
Computation of the Unweighted UniFrac diversity using only the OTU IDs present in the Tree.", log.output)
