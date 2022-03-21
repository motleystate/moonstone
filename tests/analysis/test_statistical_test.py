from unittest import TestCase

from math import isnan
import numpy as np
from packaging import version
import pandas as pd
import scipy

from moonstone.analysis.statistical_test import (
    statistical_test_groups_comparison,
    mann_whitney_u,
    ttest_independence,
    _compute_max_nbins_contingency_table,
    _compute_contingency_table,
    _comparison_chi2_pval_depending_on_nbins,
    compute_contingency_table,
    chi2_contingency
)
# from moonstone.analysis.statistical_test import *


class TestStatisticalTestFunction(TestCase):

    def setUp(self):
        self.test_df = pd.Series({
            'sample1': 0.72,
            'sample2': 1.98,
            'sample3': 0.00,
            'sample4': 0.00,
            'sample5': 4.99
        })

    def test_mann_whitney_u_groups(self):
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

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u')
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_groups_wrong_metadata(self):
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1]
            ],
            columns=['sex', 'group'],
            index=['sample6', 'sample7'],
        )

        with self.assertRaises(RuntimeError) as cm:
            statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u')
        the_exception = cm.exception
        self.assertEqual(the_exception.__str__(), "All groups have been dropped: not enough observations by group.")

    def test_mann_whitney_u_groups_nonsym(self):
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
                [np.nan, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )
        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u',
                                                    output='dataframe', sym=False)
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_groups_series_output(self):
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

        expected_df = pd.Series({(1, 2): 0.3742593192802244, (2, 1): 0.3742593192802244})
        expected_df.index.names = ['Group1', 'Group2']

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u',
                                                    output='series', sym=True)
        pd.testing.assert_series_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_groups_nonsym_series_output(self):
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

        expected_df = pd.Series({(1, 2): 0.3742593192802244})
        expected_df.index.names = ['Group1', 'Group2']

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='mann_whitney_u',
                                                    output='series', sym=False)
        pd.testing.assert_series_equal(matrix, expected_df, check_dtype=False)

    def test_mann_whitney_u_authorized_keys(self):
        """
        depending on the version of scipy check different things
        """
        with self.assertLogs("moonstone.analysis.statistical_test", level="WARNING") as log:
            stat, pval = mann_whitney_u(
                self.test_df.iloc[:2], self.test_df.iloc[2:], nan_policy="propagate"
            )
            scipy_version = version.parse(scipy.__version__)
            if scipy_version < version.parse("1.8.0"):
                self.assertEqual(len(log.output), 1)
                self.assertIn(
                    f"WARNING:moonstone.analysis.statistical_test:['nan_policy'] not available with version \
{scipy_version} of scipy.",
                    log.output,
                )
            self.assertEqual(stat, 4.0)
            self.assertEqual(pval, 0.7670968684102772)

    def test_ttest_independance_groups(self):
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
                [np.nan, 0.483379],
                [0.483379, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )

        matrix = statistical_test_groups_comparison(self.test_df, metadata_df['group'], stat_test='ttest_independence')
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)

    def test_ttest_independence(self):
        """
        depending on the version of scipy check different things
        """
        with self.assertLogs("moonstone.analysis.statistical_test", level="WARNING") as log:
            stat, pval = ttest_independence(
                self.test_df.iloc[:2], self.test_df.iloc[2:], trim=0
            )
            scipy_version = version.parse(scipy.__version__)
            if scipy_version < version.parse("1.7.0"):
                self.assertEqual(len(log.output), 1)
                self.assertIn(
                    f"WARNING:moonstone.analysis.statistical_test:['trim'] not available with version {scipy_version} \
of scipy.",
                    log.output,
                )
            self.assertEqual(stat, -0.176164072157316)
            self.assertEqual(pval, 0.8733823222145742)


class TestChi2Functions(TestCase):

    def setUp(self):
        self.test_df = pd.Series([7, 11, 13, 18, 14, 17, 9, 0, 8, 17, 8, 19, 8, 11, 1, 6, 16, 16, 15, 18, 6, 4, 3, 4,
                                  18, 16, 8, 2, 7, 14, 11, 9, 15, 10, 9, 9, 1, 7, 14, 18, 14, 12, 18, 14, 8, 4, 3, 7,
                                  2, 5, 6, 3, 19, 19, 8, 12, 7, 1, 8, 1, 11, 14, 17, 13, 10, 6, 11, 1, 5, 18, 1, 11, 1,
                                  16, 17, 19, 9, 0, 17, 7, 18, 0, 15, 7, 2, 12, 7, 9, 1, 6, 15, 9, 5, 19, 13, 17, 1, 17,
                                  3, 0])
        self.metadata_df = pd.Series([3, 2, 1, 1, 1, 2, 3, 2, 2, 2, 2, 3, 1, 3, 1, 3, 2, 3, 1, 1, 3, 1, 2, 1, 1, 1, 1,
                                      3, 3, 2, 2, 2, 2, 2, 3, 2, 3, 3, 1, 2, 1, 3, 1, 2, 1, 1, 1, 1, 3, 3, 2, 2, 3, 3,
                                      2, 2, 2, 1, 1, 1, 1, 3, 1, 1, 3, 2, 2, 3, 1, 2, 1, 2, 1, 2, 3, 2, 2, 3, 3, 1, 3,
                                      3, 3, 1, 3, 2, 2, 2, 1, 2, 3, 2, 1, 1, 3, 2, 3, 1, 1, 2])

    def test__compute_contingency_table_binsedges(self):
        # testing cut_type == "equal-width" + force_computation + warning
        expected_tab = pd.DataFrame(
            [
                [20, 17, 16],
                [15, 18, 14]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.001, 9.0), (9.0, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        expected_tab.index = ci
        expected_tab.columns.name = "col_0"
        tab = _compute_contingency_table(
            self.test_df,
            self.metadata_df,
            bins=[-0.001, 9.0, 19.0],
            cut_type="equal-width",
        )
        pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)

    def test__compute_contingency_table_binsedges_not_including_every_value(self):
        expected_tab = pd.DataFrame(
            [
                [14, 15, 8],
                [10, 15, 10]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(2, 9.0), (9.0, 17.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        expected_tab.index = ci
        expected_tab.columns.name = "col_0"
        tab = _compute_contingency_table(
            self.test_df,
            self.metadata_df,
            bins=[2, 9.0, 17.0],
            cut_type="equal-width",
        )
        pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)

    def test__compute_contingency_table_binsedges_lessthan5(self):
        # testing a cell with less than 5 observations
        # if force_coputation == True,
        # expected_tab = pd.DataFrame([
        #    [6, 2, 5],       # (-0.001, 1.0]
        #    [14, 15, 11],    # (1.0, 9.0]
        #    [15, 18, 14]     # (9.0, 19.0]
        # ], columns=[1, 2, 3])
        with self.assertLogs("moonstone.analysis.statistical_test", level="WARNING") as log:
            tab = _compute_contingency_table(
                self.test_df,
                self.metadata_df,
                bins=[-0.001, 1.0, 9.0, 19.0],
                cut_type="equal-width",
            )
            self.assertTrue(isnan(tab))
            self.assertEqual(len(log.output), 1)
            self.assertIn(
                "WARNING:moonstone.analysis.statistical_test:Some cells have less than 5 observations. \
Another statistical test would be more appropriate to compare these groups.",
                log.output,
            )

    def test__compute_contingency_table_equalwidth(self):
        # testing cut_type="equal-width" + force_computation + warning
        expected_tab = pd.DataFrame(
            [
                [11, 4, 8],
                [9, 13, 8],
                [6, 9, 5],
                [9, 9, 9]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 4.75), (4.75, 9.5), (9.5, 14.25), (14.25, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        expected_tab.index = ci
        expected_tab.columns.name = "col_0"
        with self.assertLogs("moonstone.analysis.statistical_test", level="WARNING") as log:
            tab = _compute_contingency_table(
                self.test_df,
                self.metadata_df,
                bins=4,
                cut_type="equal-width",
                force_computation=True
            )
            pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)
            self.assertEqual(len(log.output), 1)
            self.assertIn(
                "WARNING:moonstone.analysis.statistical_test:Some cells have less than 5 observations. \
Another statistical test would be more appropriate to compare these groups.",
                log.output,
            )

    def test__compute_contingency_table_equalsize(self):
        # testing cut_type="equal-size" + force_computation
        expected_tab = pd.DataFrame(
            [
                [13, 4, 9],
                [7, 13, 7],
                [7, 10, 7],
                [8, 8, 7]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.001, 5.0), (5.0, 9.0), (9.0, 15.0), (15.0, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        expected_tab.index = ci
        expected_tab.columns.name = "col_0"
        tab = _compute_contingency_table(
            self.test_df,
            self.metadata_df,
            bins=4,
            cut_type="equal-size",
            force_computation=True,
            warn=False
        )
        pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)

    def test_compute_max_nbins_contingency_table(self):
        expected_tab = pd.DataFrame(
            [
                [13, 7, 11],
                [8, 17, 8],
                [14, 11, 11]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 6.333), (6.333, 12.667), (12.667, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        expected_tab.index = ci
        expected_tab.columns.name = "col_0"
        tab, nbins = _compute_max_nbins_contingency_table(
            self.test_df,
            self.metadata_df,
        )
        pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)
        self.assertEqual(nbins, 3)

    def test_comparison_chi2_pval_depending_on_nbins(self):
        tab_nbins3 = pd.DataFrame(
            [
                [13, 7, 11],
                [8, 17, 8],
                [14, 11, 11]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 6.333), (6.333, 12.667), (12.667, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        tab_nbins3.index = ci
        tab_nbins3.columns.name = "col_0"
        # -------
        tab_nbins2 = pd.DataFrame(
            [
                [20, 17, 16],
                [15, 18, 14]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 9.5), (9.5, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "row_0"
        tab_nbins2.index = ci
        tab_nbins2.columns.name = "col_0"
        # -------
        dic_df = {
            3: [
                6.492769, 0.165247, 4,
                np.array([
                    [10.85, 10.85,  9.3],
                    [11.55, 11.55,  9.9],
                    [12.6, 12.6, 10.8]]),
                tab_nbins3
            ],
            2: [
                0.518055, 0.771802, 2,
                np.array([
                    [18.55, 18.55, 15.9],
                    [16.45, 16.45, 14.1]]),
                tab_nbins2
            ]
        }
        expected_comparison_df = pd.DataFrame.from_dict(
            dic_df, orient='index',
            columns=['chi2', 'pval', 'dof', 'expected table', 'observed table']
            )
        # -------
        comparison_df = _comparison_chi2_pval_depending_on_nbins(
            self.test_df,
            self.metadata_df,
        )
        pd.testing.assert_frame_equal(comparison_df, expected_comparison_df, check_dtype=False)

    def test_compute_contingency_table_max_nbins_na_force_computation(self):
        # testing na=True when na in numerical_series
        # + warning in _compute_max_nbins_contingency_table when force_computation=True
        expected_tab = pd.DataFrame(
            [
                [20, 16, 15],
                [15, 17, 13],
                [0, 2, 2]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(0.982, 10.0), (10.0, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci = ci.add_categories("NaN")
        ci = ci.append(pd.Index(["NaN"]))
        ci.name = "numerical_series"
        expected_tab.index = ci
        expected_tab.columns.name = "categorical_series"
        with self.assertLogs("moonstone.analysis.statistical_test", level="WARNING") as log:
            tab = compute_contingency_table(
                self.test_df.replace(0, np.nan),
                self.metadata_df,
                bins="max nbins",
                na=True,
                force_computation=True
            )
            pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)
            self.assertEqual(len(log.output), 1)
            self.assertIn(
                "WARNING:moonstone.analysis.statistical_test:moonstone wasn't able to compute a contingency table \
with at least 5 observations per cell. Another statistical test would be more appropriate to compare these groups.\n\
force_computation return contingency table of 2 x 3.",
                log.output,
            )

    def test_compute_contingency_table_max_nbins_na_not_force_computation(self):
        # testing na=True when na in categorical_series
        # + warning in _compute_max_nbins_contingency_table when force_computation=False
        with self.assertLogs("moonstone.analysis.statistical_test", level="WARNING") as log:
            tab = compute_contingency_table(
                self.test_df,
                pd.Series(
                    np.where(self.test_df == 0, np.nan, self.metadata_df),
                    index=self.metadata_df.index),
                bins="max nbins",
                na=True,
                force_computation=False
            )
            self.assertTrue(isnan(tab))
            self.assertEqual(len(log.output), 1)
            self.assertIn(
                "WARNING:moonstone.analysis.statistical_test:moonstone wasn't able to compute a contingency table \
with at least 5 observations per cell. Another statistical test would be more appropriate to compare these groups.",
                log.output,
            )

    def test_compute_contingency_table_best_pvalue(self):
        expected_tab = pd.DataFrame(
            [
                [13, 7, 11],
                [8, 17, 8],
                [14, 11, 11]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 6.333), (6.333, 12.667), (12.667, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "numerical_series"
        expected_tab.index = ci
        expected_tab.columns.name = "categorical_series"
        tab = compute_contingency_table(
            self.test_df,
            self.metadata_df,
            bins="best pvalue"
        )
        pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)

    def test_compute_contingency_given_nbins(self):
        expected_tab = pd.DataFrame(
            [
                [13, 7, 11],
                [8, 17, 8],
                [14, 11, 11]
            ],
            columns=[1, 2, 3],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 6.333), (6.333, 12.667), (12.667, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "numerical_series"
        expected_tab.index = ci
        expected_tab.columns.name = "categorical_series"
        tab = compute_contingency_table(
            self.test_df,
            self.metadata_df,
            bins=3
        )
        pd.testing.assert_frame_equal(tab, expected_tab, check_dtype=False)

    def test_chi2_contingency_best_pvalue(self):
        expected_observed_tab = pd.DataFrame(
            [
                [13, 7],
                [8, 17],
                [14, 11]
            ],
            columns=['series1', 'series2'],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 6.333), (6.333, 12.667), (12.667, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "number"
        expected_observed_tab.index = ci
        expected_observed_tab.columns.name = "category"
        chi2, pval, dof, expected_by_chi2_table, observed_table = chi2_contingency(
            self.test_df[self.metadata_df[self.metadata_df == 1].index],
            self.test_df[self.metadata_df[self.metadata_df == 2].index],
            bins="best pvalue",
            rettab=True,
        )
        self.assertEqual(chi2, 5.3999999999999995)
        self.assertEqual(pval, 0.06720551273974977)
        self.assertEqual(dof, 2)
        np.testing.assert_array_equal(
            expected_by_chi2_table,
            np.array([[10.0, 10.0], [12.5, 12.5], [12.5, 12.5]])
            )
        pd.testing.assert_frame_equal(observed_table, expected_observed_tab, check_dtype=False)

    def test_chi2_contingency_else(self):
        """
        testing the else path with bins="max nbins"
        """
        expected_observed_tab = pd.DataFrame(
            [
                [13, 7],
                [8, 17],
                [14, 11]
            ],
            columns=['series1', 'series2'],
        )
        ii = pd.IntervalIndex.from_tuples([(-0.019, 6.333), (6.333, 12.667), (12.667, 19.0)])
        ci = pd.CategoricalIndex(ii, ordered=True)
        ci.name = "number"
        expected_observed_tab.index = ci
        expected_observed_tab.columns.name = "category"
        chi2, pval, dof, expected_by_chi2_table, observed_table = chi2_contingency(
            self.test_df[self.metadata_df[self.metadata_df == 1].index],
            self.test_df[self.metadata_df[self.metadata_df == 2].index],
            bins="max nbins",
            rettab=True,
        )
        self.assertEqual(chi2, 5.3999999999999995)
        self.assertEqual(pval, 0.06720551273974977)
        self.assertEqual(dof, 2)
        np.testing.assert_array_equal(
            expected_by_chi2_table,
            np.array([[10.0, 10.0], [12.5, 12.5], [12.5, 12.5]])
            )
        pd.testing.assert_frame_equal(observed_table, expected_observed_tab, check_dtype=False)

    def test_chi2_contingency_lessthan10_and_not_force_computation(self):
        """
        testing one of the 2 series with less than 10 observations
        (not enough to bin and have at least 5 observations per cell)
        """
        chi2, pval, dof, expected_by_chi2_table, observed_table = chi2_contingency(
            self.test_df.iloc[:5],
            self.test_df.loc[5:],
            bins="max nbins",
            rettab=True,
        )
        self.assertTrue(isnan(chi2))
        self.assertTrue(isnan(pval))
        self.assertTrue(isnan(dof))
        self.assertTrue(isnan(expected_by_chi2_table))
        self.assertTrue(isnan(observed_table))

    def test_chi2_contingency_groups(self):
        expected_df = pd.DataFrame(
            [
                [np.nan, 0.06720551273974977, 0.9310443064194389],
                [0.06720551273974977, np.nan, 0.152090],
                [0.9310443064194389, 0.152090, np.nan]
            ],
            columns=[1, 2, 3],
            index=[1, 2, 3]
        )
        matrix = statistical_test_groups_comparison(
            self.test_df, self.metadata_df, stat_test='chi2_contingency'
            )
        pd.testing.assert_frame_equal(matrix, expected_df, check_dtype=False)
