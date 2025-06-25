from unittest import TestCase

from io import StringIO
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.graph_objs._cone import Cone
from skbio import TreeNode

from moonstone.analysis.diversity.beta import (
    BrayCurtis, Jaccard, WeightedUniFrac, UnweightedUniFrac
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
        pd.testing.assert_frame_equal(
            self.tested_object_instance.pcoa.samples, expected_pcoa_samples,
        )

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
        np.testing.assert_array_almost_equal(
            result_fig["data"][1]["x"], [-0.43138279,  0.56791829], decimal=8)
        np.testing.assert_array_almost_equal(
            result_fig["data"][1]["y"], [-0.06428938, -0.02690815], decimal=8)

        # PC2 vs PC3
        result_fig = self.tested_object_instance.visualize_pcoa(metadata_df, "sex", x_pc=2, y_pc=3, show=False)
        self.assertEqual(result_fig.layout["xaxis"]["title"]["text"], "PC2")
        self.assertEqual(result_fig.layout["yaxis"]["title"]["text"], "PC3")

        self.assertEqual(result_fig["data"][1]["name"], "F")
        np.testing.assert_array_almost_equal(
            result_fig["data"][1]["x"], [-0.06428938, -0.02690815], decimal=8)
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
        np.testing.assert_array_almost_equal(
            result_fig.data[0].x, [-0.24953334,  0.07681437], decimal=8)
        np.testing.assert_array_almost_equal(
            result_fig.data[0].y, [-0.26419982,  0.31054222], decimal=8)

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
        self.assertAlmostEqual(result_fig.layout.annotations[0].x, -0.4313827918950854, places=15)
        self.assertAlmostEqual(result_fig.layout.annotations[0].y, 0.09119753638724884, places=15)
        # see explanation for the use of the AlmostEqual in test_pcoa

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
        np.testing.assert_array_almost_equal(
            result_fig["data"][1]["x"], [-0.43138279,  0.56791829], decimal=8)
        np.testing.assert_array_almost_equal(
            result_fig["data"][1]["y"], [-0.06428938, -0.02690815], decimal=8)
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
