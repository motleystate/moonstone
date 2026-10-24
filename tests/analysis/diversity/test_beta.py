from unittest import TestCase

from io import StringIO
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.graph_objs._cone import Cone
from skbio import TreeNode

from moonstone.analysis.diversity.beta import (
    BrayCurtis, Jaccard, WeightedUniFrac, UnweightedUniFrac,
    permutation_p_value_se, resolve_n_jobs, too_few_permutations_warning
)


def assert_axis_almost_equal(actual, expected, decimal=8):
    """
    Compare coordinates along one PCoA axis up to a sign flip.

    Each principal coordinate is an eigenvector, and -v is as valid as v: the sign chosen depends on the
    numerical library (NumPy, BLAS) and its version. Distances between samples, and so the analysis, are the same.
    """
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    sign = -1.0 if np.dot(actual, expected) < 0 else 1.0
    np.testing.assert_array_almost_equal(sign * actual, expected, decimal=decimal)


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
        self.tested_object_instance = BrayCurtis(self.tested_object)

    def test_compute_beta_diversity_df(self):
        expected_object = pd.DataFrame.from_dict(
            {
                'sample1': [0, 0.333, 0.714],
                'sample2': [0.333, 0, 1.0],
                'sample3': [0.714, 1.0, 0],
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )
        pd.testing.assert_frame_equal(
            self.tested_object_instance.beta_diversity_df, expected_object,
            rtol=0.01    # set relative tolerance to 0.01 so that 1/3=0.333
            # check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
        )

    def test_compute_beta_diversity_series(self):
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
            self.tested_object_instance.beta_diversity_series, expected_object,
            rtol=0.01
        )
        pd.testing.assert_series_equal(
            self.tested_object_instance.diversity_indexes, expected_object,
            rtol=0.01
        )

    def test_run_statistical_test_groups_with_NaN(self):
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

        pval = self.tested_object_instance._run_statistical_test_groups(
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

        pval = self.tested_object_instance._run_statistical_test_groups(
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
        expected_object = pd.DataFrame.from_dict(
            {
                'sample2-sample3': [1.0, 'F'],
            },
            orient='index', columns=['beta_index', 'sex']
        )
        output = self.tested_object_instance._get_grouped_df_series(metadata_ser)
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
        expected_object = pd.DataFrame.from_dict(
            {
                'sample2-sample3': [1.0, 'F'],
            },
            orient='index', columns=['beta_index', 'sex']
        )
        output = self.tested_object_instance.analyse_groups(metadata_df, 'sex', show=False, show_pval=False)
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

    def test_pcoa(self):
        # So I discovered that pcoa.samples is not identical between my moonstone environment and my jupyter notebook
        # environment at the 17th decimal. Unsure why...
        # moonstone environment: python=3.10.16 =/= jupyter notebook environment: python=3.10.14
        # In both environment: scikit-bio==0.5.9, scikit-learn==1.3.1
        expected_pcoa_samples = pd.DataFrame.from_dict(
            {
                'sample1': [-0.136535, 0.091198, 0.0],
                'sample2': [-0.431383, -0.064289, 0.0],
                'sample3': [0.567918, -0.026908, 0.0],
            },
            orient='index', columns=['PC1', 'PC2', 'PC3']
        )
        samples = self.tested_object_instance.pcoa.samples
        self.assertListEqual(list(samples.columns), list(expected_pcoa_samples.columns))
        self.assertListEqual(list(samples.index), list(expected_pcoa_samples.index))
        for axis in expected_pcoa_samples.columns:
            assert_axis_almost_equal(samples[axis], expected_pcoa_samples[axis], decimal=5)

        expected_pcoa_proportions_explained = pd.Series(
            [0.975623, 0.024377, 0], index=["PC1", "PC2", "PC3"],
        )
        pd.testing.assert_series_equal(
            self.tested_object_instance.pcoa.proportion_explained, expected_pcoa_proportions_explained,
        )

    def test_label_with_proportion(self):
        expected_object = "PC1 (97.56%)"
        self.assertEqual(self.tested_object_instance._label_with_proportion("PC1"), expected_object)

    def test_unit_scale_with_scale_between0and1(self):
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
            self.tested_object_instance._unit_scale(tested_object, 228, 0.456),
            expected_object
        )

    def test_unit_scale_with_scale_above1(self):
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
            self.tested_object_instance._unit_scale(tested_object, 228, 22.8),
            expected_object
        )

    def test_scale_biplot(self):
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
            self.tested_object_instance._scale_biplot(tested_features_df, tested_samples_df),
            expected_object
        )

    def test_map_loadings(self):
        # already sorted
        tested_pci = [-18, -12, -0.8, 0, 1.3, 4.5, 12, 38.2]
        expected_object = ([6, 7], [0, 1])
        self.assertTupleEqual(self.tested_object_instance._map_loadings(tested_pci, 2), expected_object)

        # unsorted
        tested_pci = [1.3, -0.8, 38.2, 0, 4.5, -18, -12, 12]
        expected_object = ([7, 2], [5, 6])
        self.assertTupleEqual(self.tested_object_instance._map_loadings(tested_pci, 2), expected_object)

    def test_bi_plotly_multiindex(self):
        # Testing if it works with MultiIndex for the 2d version
        fig = go.Figure()
        tested_pcx = [-0.44, 0.52, 0.22, 0.07]
        tested_pcy = [-0.26, -0.07, -0.04, 0.31]
        tested_features = pd.MultiIndex.from_tuples(
            [('genus1', 'species1'), ('genus1', 'species2'), ('genus2', 'species3'), ('genus2', 'species4')],
            names=['genus', 'species'])

        result_fig = self.tested_object_instance._bi_plotly(fig, tested_pcx, tested_pcy, tested_features, 1)
        # pcx explained most by:
        #   * 'species1' in the negative direction
        #   * 'species2' in the positive direction
        # pcy explained most by:
        #   * 'species1' in the negative direction
        #   * 'species4' in the positive direction
        # So 3 features shown * 2 (arrow + text) = 6 annotations
        self.assertEqual(len(result_fig.layout.annotations), 6)

        # checking for one feature
        # 0 = arrow for first feature
        self.assertTrue(result_fig.layout.annotations[0].showarrow)
        self.assertEqual(result_fig.layout.annotations[0].x, -0.44)
        self.assertEqual(result_fig.layout.annotations[0].y, -0.26)
        # 1 = text for first feature
        self.assertEqual(result_fig.layout.annotations[1].text, 'species1')   # the real test here

        # So let's check for the two other features
        self.assertEqual(result_fig.layout.annotations[3].text, 'species2')
        self.assertEqual(result_fig.layout.annotations[5].text, 'species4')

    def test_bi_plotly3d_multiindex(self):
        # Testing if it works with MultiIndex for the 3d version
        fig = go.Figure()
        tested_pcx = [-0.44, 0.52, 0.22, 0.07]
        tested_pcy = [-0.26, -0.07, -0.04, 0.31]
        tested_pcz = [-0.08, -0.02, 0.18, -0.10]
        tested_features = pd.MultiIndex.from_tuples(
            [('genus1', 'species1'), ('genus1', 'species2'), ('genus2', 'species3'), ('genus2', 'species4')],
            names=['genus', 'species'])

        result_fig = self.tested_object_instance._bi_plotly3d(
            fig, tested_pcx, tested_pcy, tested_pcz, tested_features, 1)
        # pcx explained most by:
        #   * 'species1' in the negative direction
        #   * 'species2' in the positive direction
        # pcy explained most by:
        #   * 'species1' in the negative direction
        #   * 'species4' in the positive direction
        # pcz explained most by:
        #   * 'species4' in the negative direction
        #   * 'species3' in the positive direction
        # So 4 features shown * 3 (line + cone + text) = 12 annotations
        self.assertEqual(len(result_fig.data), 12)
        # checking for one feature
        # 0 = line for first feature
        self.assertTrue(result_fig.data[0].mode, 'lines')
        self.assertTrue(result_fig.data[0].x, [0, -0.44])
        self.assertTrue(result_fig.data[0].y, [0, -0.26])
        self.assertTrue(result_fig.data[0].z, [0, -0.08])
        # 1 = text for first feature
        self.assertEqual(result_fig.data[1].text[0], 'species1')
        # 2 = cone for first feature
        self.assertEqual(type(result_fig.data[2]), Cone)

        # Let's check for the three other features' name
        self.assertEqual(result_fig.data[4].text[0], 'species2')
        self.assertEqual(result_fig.data[7].text[0], 'species3')
        self.assertEqual(result_fig.data[10].text[0], 'species4')

    def test_visualize_pcoa_scatter_1groupcol_nobiplot_noproportion(self):
        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M'],
                'sample2': ['F'],
                'sample3': ['F'],
            },
            orient='index', columns=['sex']
        )

        # PC1 vs PC2 (default)
        # testing that it runs without error and checking a few variables
        result_fig = self.tested_object_instance.visualize_pcoa(metadata_df, "sex", show=False)
        self.assertEqual(result_fig.data[0]["type"], "scatter")
        self.assertEqual(result_fig.layout["xaxis"]["title"]["text"], "PC1")
        self.assertEqual(result_fig.layout["yaxis"]["title"]["text"], "PC2")

        self.assertEqual(result_fig["data"][1]["name"], "F")
        assert_axis_almost_equal(result_fig["data"][1]["x"], [-0.43138279,  0.56791829], decimal=8)
        assert_axis_almost_equal(result_fig["data"][1]["y"], [-0.06428938, -0.02690815], decimal=8)

        # PC2 vs PC3
        result_fig = self.tested_object_instance.visualize_pcoa(metadata_df, "sex", x_pc=2, y_pc=3, show=False)
        self.assertEqual(result_fig.layout["xaxis"]["title"]["text"], "PC2")
        self.assertEqual(result_fig.layout["yaxis"]["title"]["text"], "PC3")

        self.assertEqual(result_fig["data"][1]["name"], "F")
        assert_axis_almost_equal(result_fig["data"][1]["x"], [-0.06428938, -0.02690815], decimal=8)
        np.testing.assert_array_equal(
            result_fig["data"][1]["y"], [0, 0])

    def test_visualize_pcoa_invalidmode_proportion(self):
        # Checking that proportion work but also giving an invalid mode
        # that will be changed to the default meaning 'scatter'
        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M'],
                'sample2': ['F'],
                'sample3': ['F'],
            },
            orient='index', columns=['sex']
        )

        with self.assertLogs("moonstone.analysis.diversity.beta", level="WARNING") as log:
            result_fig = self.tested_object_instance.visualize_pcoa(
                metadata_df, "sex", x_pc=1, y_pc=2,
                mode="invalidmode", proportions=True,
                show=False
                )

            # check warning
            self.assertEqual(len(log.output), 1)
            self.assertIn(
                "WARNING:moonstone.analysis.diversity.beta:'invalidmode' not a available mode, \
set to default (scatter).",
                log.output,
            )

            # check invalid mode has been changed to 'scatter' the default
            self.assertEqual(result_fig.data[0]["type"], "scatter")

            # check proportions in title
            self.assertEqual(result_fig.layout["xaxis"]["title"]["text"], "PC1 (97.56%)")
            self.assertEqual(result_fig.layout["yaxis"]["title"]["text"], "PC2 (2.44%)")

    def test_visualize_pcoa_scatter_2groupcol(self):
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

        # testing that it runs without error and checking a few variables
        result_fig = tested_object_instance.visualize_pcoa(
            metadata_df, group_col='sex', group_col2='Group', x_pc=1, y_pc=2,
            show=False)

        # 2 "sex" * 2 "group" = 4 categories
        self.assertEqual(len(result_fig.data), 4)
        self.assertEqual(result_fig.data[0].marker.symbol, 0)
        self.assertEqual(result_fig.data[0].marker.color, '#A63A50')
        self.assertEqual(result_fig.data[0].name, "M - A")
        assert_axis_almost_equal(result_fig.data[0].x, [-0.24953334,  0.07681437], decimal=8)
        assert_axis_almost_equal(result_fig.data[0].y, [-0.26419982,  0.31054222], decimal=8)

    def test_visualize_pcoa_scatter_biplot(self):
        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M'],
                'sample2': ['F'],
                'sample3': ['F'],
            },
            orient='index', columns=['sex']
        )

        # testing that it runs without error and checking a few variables
        result_fig = self.tested_object_instance.visualize_pcoa(
            metadata_df, "sex", x_pc=1, y_pc=2,
            n_biplot_features=1,
            show=False)

        # x = pc1 axis explained most by:
        #   * 'species1' in the negative direction
        #   * 'species2' in the positive direction
        # y = pc2 axis explained most by:
        #   * 'species4' in the negative direction
        #   * 'species1' in the positive direction
        # So 3 features shown * 2 (arrow + text) = 6 annotations
        self.assertEqual(len(result_fig["layout"]["annotations"]), 6)

        # checking one annotation
        self.assertTrue(result_fig.layout.annotations[0].showarrow)
        # arrows follow the orientation of their axis, so compare up to a sign flip (see assert_axis_almost_equal)
        self.assertAlmostEqual(abs(result_fig.layout.annotations[0].x), 0.4313827918950854, places=15)
        self.assertAlmostEqual(abs(result_fig.layout.annotations[0].y), 0.09119753638724884, places=15)

    def test_visualize_pcoa_scatter3d_1groupcol_nobiplot_noproportion(self):
        metadata_df = pd.DataFrame.from_dict(
            {
                'sample1': ['M'],
                'sample2': ['F'],
                'sample3': ['F'],
            },
            orient='index', columns=['sex']
        )

        # PC1 vs PC2 vs PC3 (default)
        # testing that it runs without error and checking a few variables
        result_fig = self.tested_object_instance.visualize_pcoa(metadata_df, "sex", show=False, mode="scatter3d")
        self.assertEqual(result_fig.data[0]["type"], "scatter3d")
        self.assertEqual(result_fig.layout["scene"]["xaxis"]["title"]["text"], "PC1")
        self.assertEqual(result_fig.layout["scene"]["yaxis"]["title"]["text"], "PC2")
        self.assertEqual(result_fig.layout["scene"]["zaxis"]["title"]["text"], "PC3")

        self.assertEqual(result_fig["data"][1]["name"], "F")
        assert_axis_almost_equal(result_fig["data"][1]["x"], [-0.43138279,  0.56791829], decimal=8)
        assert_axis_almost_equal(result_fig["data"][1]["y"], [-0.06428938, -0.02690815], decimal=8)
        np.testing.assert_array_equal(
            result_fig["data"][1]["z"], [0.0, 0.0])

        # PC3 vs PC1 vs PC2
        result_fig = self.tested_object_instance.visualize_pcoa(
            metadata_df, "sex", show=False, mode="scatter3d",
            x_pc=3, y_pc=1, z_pc=2
        )
        self.assertEqual(result_fig.layout["scene"]["xaxis"]["title"]["text"], "PC3")
        self.assertEqual(result_fig.layout["scene"]["yaxis"]["title"]["text"], "PC1")
        self.assertEqual(result_fig.layout["scene"]["zaxis"]["title"]["text"], "PC2")

        np.testing.assert_array_equal(
            result_fig["data"][1]["x"], [0.0, 0.0])

    def test_visualize_pcoa_scatter3d_2groupcol_biplot_proportion(self):
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

        # testing that it runs without error and checking a few variables
        result_fig = tested_object_instance.visualize_pcoa(
            metadata_df, group_col='sex', group_col2="Group",
            mode='scatter3d', n_biplot_features=1, proportions=True,
            show=False)

        # x = pc1 axis explained most by:
        #   * 'species1' in the negative direction
        #   * 'species2' in the positive direction
        # y = pc2 axis explained most by:
        #   * 'species1' in the negative direction
        #   * 'species4' in the positive direction
        # z = pc3 axis explained most by:
        #   * 'species4' in the negative direction
        #   * 'species3' in the positive direction
        # So 4 features shown * 3 (line + cone + text) + 4 groups ("M - A", "M - B", "F - A", "F - B") = 16
        self.assertEqual(len(result_fig.data), 16)
        self.assertEqual(result_fig.layout.scene.xaxis.title.text, 'PC1 (67.37%)')
        self.assertEqual(result_fig.layout.scene.yaxis.title.text, 'PC2 (26.23%)')
        self.assertEqual(result_fig.layout.scene.zaxis.title.text, 'PC3 (5.59%)')


class TestJaccard(TestCase):

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
        tested_object_instance = Jaccard(self.tested_object)
        expected_object = pd.DataFrame.from_dict(
            {
                'sample1': [0, 0.666667, 0.5],
                'sample2': [0.666667, 0, 1],
                'sample3': [0.5, 1, 0],
            },
            orient='index', columns=['sample1', 'sample2', 'sample3']
        )
        pd.testing.assert_frame_equal(
            tested_object_instance.beta_diversity_df, expected_object,
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


class TestPermanova(TestCase):
    """PERMANOVA and PERMDISP on synthetic counts with known structure."""

    @staticmethod
    def make_counts(spread_a, spread_b, shift_b=1.0, n=15, n_taxa=20, seed=1):
        rng = np.random.default_rng(seed)
        base = rng.gamma(2, 50, size=n_taxa)
        shifted = base.copy()
        shifted[:5] *= shift_b
        a = [rng.poisson(base * rng.lognormal(0, spread_a, n_taxa)) for _ in range(n)]
        b = [rng.poisson(shifted * rng.lognormal(0, spread_b, n_taxa)) for _ in range(n)]
        samples = [f"s{i}" for i in range(2 * n)]
        counts = pd.DataFrame(np.vstack(a + b).T, index=[f"t{i}" for i in range(n_taxa)], columns=samples)
        metadata = pd.DataFrame({"group": ["A"] * n + ["B"] * n}, index=samples)
        return counts, metadata

    def test_location_difference_detected(self):
        counts, metadata = self.make_counts(0.2, 0.2, shift_b=8.0)
        res = BrayCurtis(counts).permanova(metadata, "group")
        self.assertLess(res["p_value"], 0.01)
        self.assertEqual(res["n_samples"], 30)
        self.assertEqual(res["n_groups"], 2)
        self.assertEqual(res.name, "group")

    def test_dispersion_difference_flagged_by_permdisp(self):
        counts, metadata = self.make_counts(0.1, 1.0)
        res = BrayCurtis(counts).permanova(metadata, "group")
        self.assertLess(res["permdisp_p_value"], 0.01)

    def test_no_dispersion(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        res = BrayCurtis(counts).permanova(metadata, "group", dispersion=False)
        self.assertTrue(np.isnan(res["permdisp_p_value"]))

    def test_seed_reproducible(self):
        counts, metadata = self.make_counts(0.3, 0.3)
        metadata["random"] = np.random.default_rng(3).choice(["x", "y"], 30)
        div = BrayCurtis(counts)
        first = div.permanova(metadata, "random", permutations=99, seed=5)
        second = div.permanova(metadata, "random", permutations=99, seed=5)
        pd.testing.assert_series_equal(first, second)

    def test_r2_by_hand(self):
        counts = pd.DataFrame(
            {"s1": [10, 0, 5], "s2": [8, 1, 6], "s3": [9, 2, 4], "s4": [0, 10, 3], "s5": [1, 9, 2], "s6": [2, 7, 5]},
            index=["t1", "t2", "t3"],
        )
        metadata = pd.DataFrame({"g": ["a", "a", "a", "b", "b", "b"]}, index=counts.columns)
        div = BrayCurtis(counts)
        res = div.permanova(metadata, "g", permutations=9)
        d2 = div.beta_diversity.data ** 2
        ss_total = d2[np.triu_indices(6, 1)].sum() / 6
        ss_within = sum(d2[np.ix_(idx, idx)][np.triu_indices(3, 1)].sum() / 3 for idx in ([0, 1, 2], [3, 4, 5]))
        self.assertAlmostEqual(res["R2"], 1 - ss_within / ss_total)
        self.assertAlmostEqual(res["pseudo_F"], ((ss_total - ss_within) / 1) / (ss_within / 4))

    def test_missing_and_unmatched_samples_left_out(self):
        counts, metadata = self.make_counts(0.2, 0.2, shift_b=8.0)
        metadata.loc["s0", "group"] = np.nan
        metadata = metadata.drop("s1")
        metadata.loc["not_sequenced", "group"] = "A"
        res = BrayCurtis(counts).permanova(metadata, "group", permutations=99)
        self.assertEqual(res["n_samples"], 28)

    def test_missing_values_counted_not_grouped(self):
        counts, metadata = self.make_counts(0.2, 0.2, shift_b=8.0)
        complete = [s for s in counts.columns if s not in ("s0", "s1", "s15", "s16")]
        expected = BrayCurtis(counts[complete]).permanova(metadata.loc[complete], "group", permutations=99)
        missing_values = {
            "nan": np.nan, "None": None, "pd.NA": pd.NA, "empty": "", "space": " ",
        }
        for dtype in [object, "string", "category"]:
            for name, value in missing_values.items():
                with self.subTest(dtype=dtype, missing=name):
                    meta = metadata.copy()
                    meta["group"] = meta["group"].astype(object)
                    meta.loc[["s0", "s1", "s15", "s16"], "group"] = value
                    meta["group"] = meta["group"].astype(dtype)
                    res = BrayCurtis(counts).permanova(meta, "group", permutations=99)
                    self.assertEqual(res["n_samples"], 26)
                    self.assertEqual(res["n_missing"], 4)
                    self.assertEqual(res["n_groups"], 2)
                    self.assertEqual(res["pseudo_F"], expected["pseudo_F"])
                    self.assertEqual(res["permdisp_F"], expected["permdisp_F"])

    def test_missing_values_in_numeric_column(self):
        counts, metadata = self.make_counts(0.2, 0.2, shift_b=8.0)
        for dtype, value in [(float, np.nan), ("Float64", pd.NA)]:
            with self.subTest(dtype=dtype):
                meta = pd.DataFrame({"visit": [1.0] * 15 + [2.0] * 15}, index=metadata.index).astype(dtype)
                meta.loc[["s0", "s15", "s16"], "visit"] = value
                res = BrayCurtis(counts).permanova(meta, "visit", permutations=9)
                self.assertEqual(res["n_samples"], 27)
                self.assertEqual(res["n_missing"], 3)
                self.assertEqual(res["n_groups"], 2)

    def test_missing_values_logged_once(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata.loc["s0", "group"] = np.nan
        with self.assertLogs("moonstone.analysis.diversity.beta", level="WARNING") as log:
            BrayCurtis(counts).permanova_screen(metadata, permutations=99)
        self.assertEqual(log.output, ["WARNING:moonstone.analysis.diversity.beta:group: 1 sample(s) with a missing "
                                      "value are left out."])

    def test_screen_reports_missing_and_skips_with_warning(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata["mostly_missing"] = [np.nan] * 28 + ["x", "y"]
        metadata.loc["s3", "group"] = None
        with self.assertLogs("moonstone.analysis.diversity.beta", level="WARNING") as log:
            screen = BrayCurtis(counts).permanova_screen(metadata, permutations=9)
        self.assertEqual(screen.loc["group", "n_missing"], 1)
        self.assertEqual(screen.loc["group", "n_samples"], 29)
        self.assertEqual(screen.loc["mostly_missing", "n_missing"], 28)
        self.assertEqual(screen.loc["mostly_missing", "skipped"], "fewer than 2 groups with at least 2 samples")
        self.assertTrue(any("mostly_missing: not tested" in line for line in log.output))

    def test_mostly_missing_variable_warns_without_filling(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata["half"] = ["x", "y"] * 7 + [np.nan] * 16
        with self.assertLogs("moonstone.analysis.diversity.beta", level="WARNING") as log:
            res = BrayCurtis(counts).permanova(metadata, "half", permutations=9)
        self.assertTrue(any("16 of 30 samples (53%) have no value" in line for line in log.output))
        self.assertEqual(res["n_samples"], 14)
        self.assertEqual(res["n_missing"], 16)
        self.assertEqual(res["n_groups"], 2)

    def test_half_missing_does_not_warn_as_mostly_missing(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata["half"] = ["x", "y"] * 8 + [np.nan] * 14
        with self.assertLogs("moonstone.analysis.diversity.beta", level="WARNING") as log:
            BrayCurtis(counts).permanova(metadata, "half", permutations=9)
        self.assertFalse(any("have no value" in line for line in log.output))

    def test_unused_categories_not_reported(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata["group"] = pd.Categorical(metadata["group"], categories=["A", "B", "unused"])
        with self.assertNoLogs("moonstone.analysis.diversity.beta", level="INFO"):
            res = BrayCurtis(counts).permanova(metadata, "group", permutations=9)
        self.assertEqual(res["n_groups"], 2)

    def test_small_groups_left_out(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata.loc["s0", "group"] = "C"
        res = BrayCurtis(counts).permanova(metadata, "group", permutations=9, min_group_size=2)
        self.assertEqual(res["n_groups"], 2)
        self.assertEqual(res["n_samples"], 29)

    def test_one_group_raises(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata["group"] = "A"
        with self.assertRaises(ValueError):
            BrayCurtis(counts).permanova(metadata, "group")

    def test_screen_corrects_and_skips(self):
        from statsmodels.stats.multitest import multipletests
        counts, metadata = self.make_counts(0.2, 0.2, shift_b=8.0)
        rng = np.random.default_rng(0)
        for i in range(5):
            metadata[f"null{i}"] = rng.choice(["x", "y"], 30)
        metadata["constant"] = "same"
        metadata["continuous"] = rng.normal(size=30)
        screen = BrayCurtis(counts).permanova_screen(metadata, permutations=199)
        self.assertEqual(list(screen.index), list(metadata.columns))
        self.assertEqual(screen.loc["constant", "skipped"], "fewer than 2 groups with at least 2 samples")
        self.assertEqual(screen.loc["continuous", "skipped"], "more than 15 groups")
        tested = screen[screen["skipped"] == ""]
        self.assertEqual(len(tested), 6)
        np.testing.assert_allclose(
            tested["p_value_corrected"], multipletests(tested["p_value"], method="fdr_bh")[1]
        )
        self.assertTrue(screen.loc[["constant", "continuous"], "p_value_corrected"].isna().all())
        self.assertLess(screen.loc["group", "p_value_corrected"], 0.05)

    def test_screen_without_correction(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        screen = BrayCurtis(counts).permanova_screen(metadata, correction_method=None, permutations=9)
        self.assertNotIn("p_value_corrected", screen.columns)

    def test_screen_parallel_matches_serial(self):
        counts, metadata = self.make_counts(0.2, 0.4, shift_b=2.0)
        rng = np.random.default_rng(0)
        for i in range(3):
            metadata[f"null{i}"] = rng.choice(["x", "y", "z"], 30)
        metadata["constant"] = "same"
        div = BrayCurtis(counts)
        serial = div.permanova_screen(metadata, permutations=199, n_jobs=1)
        parallel = div.permanova_screen(metadata, permutations=199, n_jobs=3)
        pd.testing.assert_frame_equal(serial, parallel)

    def test_screen_warns_when_permutations_too_few(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        metadata["other"] = ["x", "y"] * 15
        with self.assertLogs("moonstone.analysis.diversity.beta", level="WARNING") as log:
            BrayCurtis(counts).permanova_screen(metadata, permutations=19)
        self.assertTrue(any("19 permutations are too few for 2 tests" in line for line in log.output))

    def test_default_permutations_and_se(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        res = BrayCurtis(counts).permanova(metadata, "group", dispersion=False)
        self.assertEqual(res["permutations"], 9999)
        self.assertAlmostEqual(res["p_value_se"], np.sqrt(res["p_value"] * (1 - res["p_value"]) / 9999))
        self.assertTrue(np.isnan(res["permdisp_p_value_se"]))

    def test_permutation_helpers(self):
        self.assertAlmostEqual(permutation_p_value_se(0.05, 9999), np.sqrt(0.05 * 0.95 / 9999))
        self.assertTrue(np.isnan(permutation_p_value_se(0.5, 0)))
        # 1 / (N + 1) times the number of tests must be below 0.05
        self.assertEqual(too_few_permutations_warning(499, 9999), "")
        self.assertIn("too few for 500 tests", too_few_permutations_warning(500, 9999))
        self.assertIn("too few for 50 tests", too_few_permutations_warning(50, 999))
        self.assertIn("Use at least 1000 permutations", too_few_permutations_warning(50, 999))
        self.assertEqual(too_few_permutations_warning(0, 9), "")
        self.assertEqual(resolve_n_jobs(3), 3)
        self.assertGreaterEqual(resolve_n_jobs(0), 1)
        self.assertEqual(resolve_n_jobs(-1), resolve_n_jobs(0))
        with self.assertRaises(ValueError):
            resolve_n_jobs(-2)

    def test_works_for_jaccard(self):
        counts, metadata = self.make_counts(0.2, 0.2)
        counts[counts < 100] = 0    # Jaccard uses presence/absence only, so some absences are needed
        res = Jaccard(counts).permanova(metadata, "group", permutations=9)
        self.assertEqual(res["n_samples"], 30)
