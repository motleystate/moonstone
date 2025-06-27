from unittest import TestCase

from io import StringIO
import numpy as np
from packaging import version
import pandas as pd
import scipy
from skbio import TreeNode

from moonstone.analysis.diversity.alpha import (
    ShannonIndex, SimpsonInverseIndex, Chao1Index,
    FaithsPhylogeneticDiversity
)


class TestShannonIndex(TestCase):

    def setUp(self):
        self.tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 0, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        self.metadata_df = pd.DataFrame(
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

    def test_compute_alpha_diversity(self):
        expected_object = pd.Series({
            'sample1': 0.721928,
            'sample2': 0.985228,
            'sample3': -0.000000,
            'sample4': -0.000000,
            'sample5': 1.992778
        })
        tested_object_instance = ShannonIndex(self.tested_object)
        pd.testing.assert_series_equal(tested_object_instance.compute_diversity(), expected_object)

    def test_visualize(self):
        # Not real test but make sure that visualize() runs without errors
        tested_object_instance = ShannonIndex(self.tested_object)
        tested_object_instance.visualize(show=False)

    def test_analyse_groups(self):
        # Simple case: only 1 group_col, pval not displayed on graph
        tested_object_instance = ShannonIndex(self.tested_object)
        expected_df = pd.DataFrame(
            [
                [np.nan, 0.374259],
                [0.374259, np.nan]
            ],
            columns=[1, 2],
            index=[1, 2]
        )

        output = tested_object_instance.analyse_groups(self.metadata_df, 'group', make_graph=False,
                                                       structure_pval='dataframe', sym=True)
        pd.testing.assert_frame_equal(output['pval'], expected_df, check_dtype=False)

    def test_pvalue_correction(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 3],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'],
        )
        expected_df = pd.DataFrame(
            [
                [np.nan, 1.0, 1.0],
                [1.0, np.nan, 1.0],
                [1.0, 1.0, np.nan]
            ],
            columns=[1, 2, 3],
            index=[1, 2, 3]
        )

        output = tested_object_instance.analyse_groups(metadata_df, 'group', make_graph=False,
                                                       structure_pval='dataframe', sym=True,
                                                       correction_method='fdr_bh'
                                                       )
        pd.testing.assert_frame_equal(output['pval'], expected_df, check_dtype=False)

    def test_invalid_correction_method_param(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        with self.assertLogs('moonstone.analysis.diversity.base', level='WARNING') as log:
            tested_object_instance._valid_correction_method_param("lalala")
            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.analysis.diversity.base:correction_method='lalala' not valid, \
set to default (None).", log.output)

    def test_valid_pval_param(self):
        # everything is valid
        # testing with none default value to be sure
        tested_object_instance = ShannonIndex(self.tested_object)
        pvc, pvd = tested_object_instance._valid_pval_param(
            "same group_col or group_col2 values", "same group_col values", group_col2=True
        )
        self.assertEqual(pvc, "same group_col or group_col2 values")
        self.assertEqual(pvd, "same group_col values")

    def test_invalid_pval_param_wogroup_col2(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        with self.assertLogs('moonstone.analysis.diversity.base', level='WARNING') as log:
            tested_object_instance._valid_pval_param("same group_col", None, group_col2=False)
            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.analysis.diversity.base:Without a second group column (group_col2), \
pval_to_compute='same group_col' not valid, set to default (all).", log.output)

    def test_invalid_pval_param_wgroup_col2(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        with self.assertLogs('moonstone.analysis.diversity.base', level='WARNING') as log:
            tested_object_instance._valid_pval_param("lalala", "lilili", group_col2=True)
            self.assertEqual(len(log.output), 2)
            self.assertIn("WARNING:moonstone.analysis.diversity.base:pval_to_compute='lalala' not valid, \
set to default (all).", log.output)
            self.assertIn("WARNING:moonstone.analysis.diversity.base:pval_to_display='lilili' not valid, \
set to default (None).", log.output)

    def test_inconsistent_pval_to_diplay_param(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        with self.assertRaises(ValueError) as cm:
            tested_object_instance._valid_pval_param("same group_col or group_col2 values", "all", group_col2=True)
        the_exception = cm.exception
        expected_msg = "pval_to_display='all' not valid, when pval_to_compute='same group_col or group_col2 values'. \
pval_to_display should be set to: ['same group_col or group_col2 values', 'same group_col values', None]"
        self.assertEqual(the_exception.__str__(), expected_msg)

    def test_analyse_groups_pval_to_compute_all(self):
        tested_object_instance = ShannonIndex(self.tested_object)

        if version.parse(scipy.__version__) > version.parse("1.6.3"):
            expected_ser = pd.Series({
                ('1 - F', '1 - M'): 1.0,
                ('1 - F', '2 - F'): 1.0,
                ('1 - F', '2 - M'): 1.0,
                ('1 - M', '2 - F'): 1.0,
                ('1 - M', '2 - M'): 0.6666666666666666,
                ('2 - F', '2 - M'): 1.0
                },
            )
        else:
            expected_ser = pd.Series({
                ('1 - F', '1 - M'): 1.0,
                ('1 - F', '2 - F'): 1.0,
                ('1 - F', '2 - M'): 1.0,
                ('1 - M', '2 - F'): 0.5402913746074199,
                ('1 - M', '2 - M'): 0.5402913746074199,
                ('2 - F', '2 - M'): 1.0
                },
            )
        # if scipy 1.6.3 = no `method` argument -> uses "asymptotic" method
        # if scipy > 1.6.3 = `method` argument; default being 'auto',
        # which chooses 'exact' when the size of one of the samples is less than 8 and there are no ties;
        # chooses 'asymptotic' otherwise.
        # -> so here, uses "exact"

        expected_ser.index.names = ["Group1", "Group2"]

        output = tested_object_instance.analyse_groups(
            self.metadata_df, group_col='group', group_col2='sex',
            make_graph=False, structure_pval='series', sym=False,
            correction_method="uncorrected",
            pval_to_compute="all"
        )
        pd.testing.assert_series_equal(output['pval'], expected_ser, check_dtype=False)

    def test_analyse_groups_pval_to_compute_same_group_col_values(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        expected_ser = pd.Series({
            ('2 - F', '2 - M'): 1.0,
            ('1 - F', '1 - M'): 1.0
            },
        )
        expected_ser.index.names = ["Group1", "Group2"]

        output = tested_object_instance.analyse_groups(
            self.metadata_df, group_col='group', group_col2='sex',
            make_graph=False, structure_pval='series', sym=False,
            pval_to_compute="same group_col values"
        )
        pd.testing.assert_series_equal(output['pval'], expected_ser, check_dtype=False)

    def test_analyse_groups_pval_to_compute_same_group_col_or_group_col2_values(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        if version.parse(scipy.__version__) > version.parse("1.6.3"):
            expected_ser = pd.Series({
                ('2 - F', '2 - M'): 1.0,
                ('1 - F', '1 - M'): 1.0,
                ('1 - F', '2 - F'): 1.0,
                ('1 - M', '2 - M'): 0.666666666666666
                },
            )
        else:
            expected_ser = pd.Series({
                ('2 - F', '2 - M'): 1.0,
                ('1 - F', '1 - M'): 1.0,
                ('1 - F', '2 - F'): 1.0,
                ('1 - M', '2 - M'): 0.5402913746074199
                },
            )
        # if scipy 1.6.3 = no `method` argument -> uses "asymptotic" method
        # if scipy > 1.6.3 = `method` argument; default being 'auto',
        # which chooses 'exact' when the size of one of the samples is less than 8 and there are no ties;
        # chooses 'asymptotic' otherwise.
        # -> so here, uses "exact"
        expected_ser.index.names = ["Group1", "Group2"]

        output = tested_object_instance.analyse_groups(
            self.metadata_df, group_col='group', group_col2='sex',
            make_graph=False, structure_pval='series', sym=False,
            pval_to_compute="same group_col or group_col2 values"
        )
        pd.testing.assert_series_equal(output['pval'], expected_ser, check_dtype=False)

    def test_pval_selection(self):
        tested_object = pd.Series({
            ('A', 'B'): 0.03,
            ('A', 'C'): 0.5,
            ('A', 'D'): 0.0014,
            ('B', 'C'): 0.2,
            ('B', 'D'): 0.001,
            ('C', 'D'): 0.00067,
        })
        tested_object.index.names = ["Group1", "Group2"]
        tested_object_instance = ShannonIndex(self.tested_object)

        expected_ser = pd.Series({
            ('A', 'B'): 0.03,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        pd.testing.assert_series_equal(
            tested_object_instance._pval_selection(tested_object, ['A', 'B', 'C']),
            expected_ser
        )

        expected_ser = pd.Series({
            ('A', 'D'): 0.0014,
            ('C', 'D'): 0.00067,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        pd.testing.assert_series_equal(
            tested_object_instance._pval_selection(tested_object, ['A', 'C', 'D']),
            expected_ser
        )

    def test_pval_selection_with_group_col2(self):
        # pval_to_compute = 'all'
        tested_object = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('M - A', 'F - A'): 0.5,
            ('M - C', 'F - A'): 0.0034,
            ('F - B', 'M - A'): 0.6,
            ('F - C', 'F - B'): 0.0014,
            ('F - B', 'M - B'): 0.2,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
            ('M - C', 'F - C'): 0.0003,
            ('M - C', 'F - B'): 0.0056,
            ('M - A', 'M - C'): 0.89,
            ('M - A', 'F - C'): 0.0006,
            ('M - B', 'F - C'): 0.0043,
            ('F - A', 'M - B'): 0.234,
        })
        tested_object.index.names = ["Group1", "Group2"]
        tested_object_instance = ShannonIndex(self.tested_object)

        groups = ['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B']

        # Case pval_to_display = 'all':  only pval > 0.05 removed
        pval_to_display = 'all'
        expected_ser = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('M - C', 'F - A'): 0.0034,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
            ('M - C', 'F - C'): 0.0003,
            ('M - C', 'F - B'): 0.0056,
            ('M - A', 'F - C'): 0.0006,
            ('M - B', 'F - C'): 0.0043,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'all', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

        # Case pval_to_display = 'same group_col or group_col2 values'
        pval_to_display = 'same group_col or group_col2 values'
        expected_ser = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
            ('M - C', 'F - C'): 0.0003,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'all', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

        # Case pval_to_display = 'same group_col values'
        pval_to_display = 'same group_col values'
        expected_ser = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'all', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

        # Case final_groups defined
        groups = ['M - A', 'F - A', 'M - C', 'F - C']
        expected_ser = pd.Series({
            ('M - C', 'F - A'): 0.0034,
            ('F - A', 'F - C'): 0.0031,
            ('M - C', 'F - C'): 0.0003,
            ('M - A', 'F - C'): 0.0006,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'all', 'all'
        )
        pd.testing.assert_series_equal(res, expected_ser)

    def test_pval_selection_with_group_col2_pval_to_compute_same_group_col_or_same_group_col2(self):
        # pval_to_compute = 'same group_col or group_col2 values'
        tested_object = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('M - A', 'F - A'): 0.5,
            ('F - C', 'F - B'): 0.0014,
            ('F - B', 'M - B'): 0.2,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
            ('F - C', 'M - C'): 0.0003,
            ('M - A', 'M - C'): 0.89,
        })
        tested_object.index.names = ["Group1", "Group2"]
        tested_object_instance = ShannonIndex(self.tested_object)

        groups = ['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B']

        # Case pval_to_display = 'same group_col or group_col2 values':  only pval > 0.05 removed
        pval_to_display = 'same group_col or group_col2 values'
        expected_ser = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
            ('F - C', 'M - C'): 0.0003,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'same group_col or group_col2 values', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

        # Case pval_to_display = 'same group_col values'
        pval_to_display = 'same group_col values'
        expected_ser = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'same group_col or group_col2 values', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

        # Case final_groups defined
        groups = ['M - A', 'F - A', 'M - C', 'F - C']
        expected_ser = pd.Series({
            ('F - A', 'F - C'): 0.0031,
            ('F - C', 'M - C'): 0.0003,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'same group_col or group_col2 values', 'same group_col or group_col2 values'
        )
        pd.testing.assert_series_equal(res, expected_ser)

    def test_pval_selection_with_group_col2_pval_to_compute_same_group_col(self):
        # pval_to_compute = 'same group_col values'
        tested_object = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
            ('M - A', 'M - C'): 0.89,
        })
        tested_object.index.names = ["Group1", "Group2"]
        tested_object_instance = ShannonIndex(self.tested_object)

        groups = ['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B']

        # Case pval_to_display = 'same group_col values'
        pval_to_display = 'same group_col values'
        expected_ser = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('F - C', 'F - B'): 0.0014,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'same group_col values', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

        # Case final_groups defined
        groups = ['M - A', 'F - A', 'M - C', 'F - C']
        expected_ser = pd.Series({
            ('F - A', 'F - C'): 0.0031,
        })
        expected_ser.index.names = ["Group1", "Group2"]
        res = tested_object_instance._pval_selection_with_group_col2(
            tested_object, groups, 'same group_col values', pval_to_display
        )
        pd.testing.assert_series_equal(res, expected_ser)

    def test_generate_ordered_final_groups(self):
        tested_object = pd.DataFrame.from_dict({
            'sex_Group': {'comp1': 'M - A', 'comp2': 'F - B', 'comp3': 'F - A', 'comp4': 'F - C', 'comp5': 'M - C',
                          'comp6': 'M - B', 'comp7': 'M - A'},
            'sex': {'comp1': 'M', 'comp2': 'F', 'comp3': 'F', 'comp4': 'F', 'comp5': 'M', 'comp6': 'M', 'comp7': 'M'},
            'Group': {'comp1': 'A', 'comp2': 'B', 'comp3': 'A', 'comp4': 'C', 'comp5': 'C', 'comp6': 'B', 'comp7': 'A'},
        })
        groups_Group = ["C", "A", "B"]
        groups_sex = ["M", "F"]
        tested_object_instance = ShannonIndex(self.tested_object)
        final_groups = tested_object_instance._generate_ordered_final_groups(
            tested_object, final_group_col='sex_Group', group_col='Group', group_col2='sex',
            groups=groups_Group, groups2=groups_sex
        )
        self.assertListEqual(final_groups, ['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B'])
        # now testing other way around
        final_groups = tested_object_instance._generate_ordered_final_groups(
            tested_object, final_group_col='sex_Group', group_col='sex', group_col2='Group',
            groups=groups_sex, groups2=groups_Group
        )
        self.assertListEqual(final_groups, ['M - C', 'M - A', 'M - B', 'F - C', 'F - A', 'F - B'])
        # testing if order of groups not given
        final_groups = tested_object_instance._generate_ordered_final_groups(
            tested_object, final_group_col='sex_Group', group_col='sex', group_col2='Group',
            groups=None, groups2=groups_Group
        )
        self.assertListEqual(final_groups, ['F - C', 'F - A', 'F - B', 'M - C', 'M - A', 'M - B'])
        # testing if order of groups2 not given
        final_groups = tested_object_instance._generate_ordered_final_groups(
            tested_object, final_group_col='sex_Group', group_col='sex', group_col2='Group',
            groups=groups_sex, groups2=None
        )
        self.assertListEqual(final_groups, ['M - A', 'M - B', 'M - C', 'F - A', 'F - B', 'F - C'])

    def test_order_pval_series(self):
        tested_object = pd.Series({
            ('M - C', 'M - B'): 0.03,
            ('M - A', 'F - A'): 0.5,
            ('F - C', 'F - B'): 0.0014,
            ('M - B', 'F - B'): 0.2,
            ('F - A', 'F - B'): 0.001,
            ('M - A', 'M - B'): 0.00067,
            ('F - A', 'F - C'): 0.0031,    # should be reorganized as ('F - C', 'F - A')
            ('M - C', 'F - C'): 0.0003,
            ('M - A', 'M - C'): 0.89,      # should be reorganized as ('M - C', 'M - A')
        })
        tested_object.index.names = ["Group1", "Group2"]
        groups = ['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B']
        tested_object_instance = ShannonIndex(self.tested_object)

        level0 = pd.Categorical(
            ['M - C', 'M - C', 'M - C', 'F - C', 'F - C', 'M - A', 'M - A', 'F - A', 'M - B'],
            categories=['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B'],
            ordered=True, dtype='category'
        )
        level1 = pd.Categorical(
            ['F - C', 'M - A', 'M - B', 'F - A', 'F - B', 'F - A', 'M - B', 'F - B', 'F - B'],
            categories=['M - C', 'F - C', 'M - A', 'F - A', 'M - B', 'F - B'],
            ordered=True, dtype='category'
        )
        data = [[0.0003, 0.89, 0.03, 0.0031, 0.0014, 0.5, 0.00067, 0.001, 0.2]]
        expected_ser = pd.DataFrame(
            data=data,
            columns=pd.MultiIndex.from_arrays([level0, level1])
        ).T[0]
        expected_ser.index.names = ["Group1", "Group2"]

        pd.testing.assert_series_equal(
            tested_object_instance._order_pval_series(tested_object, groups),
            expected_ser
        )

    def test_generate_shapes_annotations_lists_simple(self):
        tested_pval = pd.Series({
            (1, 2): 0.0005, (1, 3): 0.02, (2, 3): 0.0001, (3, 4): 0.01
            })
        tested_pval.index.names = ["Group1", "Group2"]

        tested_object_instance = ShannonIndex(self.tested_object)
        result_shapes, result_annot, result_nlevels = tested_object_instance._generate_shapes_annotations_lists_simple(
            tested_pval, [1, 2, 3, 4], 3)
        self.assertEqual(result_nlevels, 3)
        self.assertEqual(len(result_annot), 4)
        self.assertEqual(len(result_shapes), 12)  # 3 * 4 = 12

        # (1, 2)
        # from Group1 to Group2
        self.assertDictEqual(result_shapes[0], {'x0': 0, 'y0': 3.1, 'x1': 1, 'y1': 3.1, 'line': {'width': 0.53}})
        # right edge
        self.assertDictEqual(
            result_shapes[1], {'x0': 0, 'y0': 3.0700000000000003, 'x1': 0, 'y1': 3.1, 'line': {'width': 0.53}})
        # left edge
        self.assertDictEqual(
            result_shapes[2], {'x0': 1, 'y0': 3.0700000000000003, 'x1': 1, 'y1': 3.1, 'line': {'width': 0.53}})
        self.assertEqual(result_annot[0]['text'], "**")
        self.assertEqual(result_annot[0]['y'], 3.13)

        # (1, 3)
        self.assertEqual(result_annot[1]['text'], "*")
        self.assertEqual(result_annot[1]['y'], 3.23)

        # (2, 3)
        self.assertEqual(result_annot[2]['y'], 3.33)

        # (3, 4)
        self.assertEqual(result_annot[3]['text'], "**")
        self.assertEqual(result_annot[3]['y'], 3.13)

    def test_analyse_groups_pval_to_display(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 0, 0, 0, 4, 7, 0, 6],
                'species2': [1, 0, 2, 0, 5, 3, 0, 4],
                'species3': [9, 0, 0, 1, 4, 2, 0, 5],
                'species4': [8, 3, 0, 0, 9, 1, 0, 1],
                'species5': [0, 0, 0, 0, 4, 0, 0, 8],
                'species6': [4, 0, 0, 0, 4, 9, 0, 9],
                'species7': [6, 1, 0, 0, 7, 8, 0, 7],
                'species8': [1, 0, 0, 0, 8, 4, 1, 6]

            },
            orient='index',
            columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8'])

        metadata_df = pd.DataFrame(
            [
                ['F', 2],
                ['M', 1],
                ['F', 1],
                ['M', 1],
                ['M', 2],
                ['F', 2],
                ['F', 1],
                ['M', 2],
            ],
            columns=['sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8'],
        )
        tested_object_instance = ShannonIndex(tested_object)

        output = tested_object_instance.analyse_groups(
            metadata_df, 'group', make_graph=True, show=False, show_pval=False,
            pval_to_display="all", groups=[1, 2],
            structure_pval='dataframe', sym=True)
        self.assertEqual(output['fig'].layout.annotations[0].text, "*")
        self.assertEqual(len(output['fig'].layout.shapes), 3)

    def test_analyse_groups_pval_to_display_nothing_to_display(self):
        tested_object_instance = ShannonIndex(self.tested_object)

        output = tested_object_instance.analyse_groups(
            self.metadata_df, 'group', make_graph=True, show=False, show_pval=False,
            pval_to_display="all", groups=[1, 2],
            structure_pval='dataframe', sym=True)

        self.assertEqual(output['pval'].loc[1, 2], 0.3742593192802244)
        self.assertEqual(len(output['fig'].layout.shapes), 0)


class TestSimpsonInverseIndex(TestCase):

    def test_compute_alpha_diversity(self):

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 0, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = SimpsonInverseIndex(tested_object)
        expected_object = pd.Series({
            'sample1': 1.470588,
            'sample2': 1.960000,
            'sample3': 1.000000,
            'sample4': 1.000000,
            'sample5': 3.958904
        })
        pd.testing.assert_series_equal(tested_object_instance.compute_diversity(), expected_object)

    def test_visualize(self):
        # Not real test but make sure that visualize() runs without errors
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 0, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = SimpsonInverseIndex(tested_object)
        tested_object_instance.visualize(show=False)


class TestChao1Index(TestCase):

    def test_compute_alpha_diversity(self):

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 2, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = Chao1Index(tested_object, bias_corrected=False)
        expected_object = pd.Series({
            'sample1': 2.0,
            'sample2': 2.0,
            'sample3': 1.0,
            'sample4': 2.5,
            'sample5': 4.0
        })
        pd.testing.assert_series_equal(
            tested_object_instance.compute_diversity(), expected_object
            )

    def test_visualize(self):
        # Not real test but make sure that visualize() runs without errors
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 0, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = Chao1Index(tested_object)
        tested_object_instance.visualize(show=False)


class TestFaithsPhylogeneticDiversity(TestCase):

    def setUp(self):
        self.tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 2, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])

    def test_compute_alpha_diversity(self):
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species4:0.5,species5:0.5):1.0)root;'))
        tested_object_instance = FaithsPhylogeneticDiversity(self.tested_object, tree, validate=False)
        expected_object = pd.Series({
            'sample1': 1.75,
            'sample2': 3.00,
            'sample3': 1.50,
            'sample4': 2.50,
            'sample5': 4.25
        })
        pd.testing.assert_series_equal(
            tested_object_instance.compute_diversity(), expected_object
            )

    def test_compute_alpha_diversity_incomplete_tree(self):
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species6:0.5,species5:0.5):1.0)root;'))
        tested_object_instance = FaithsPhylogeneticDiversity(self.tested_object, tree, validate=False)

        with self.assertRaises(RuntimeError) as cm:
            tested_object_instance.compute_diversity()
        the_exception = cm.exception
        self.assertEqual(the_exception.__str__(), "INCOMPLETE TREE: missing ['species4'].")

    def test_compute_alpha_diversity_force_computation(self):
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
        tested_object_instance = FaithsPhylogeneticDiversity(tested_object, tree, force_computation=True)
        expected_object = pd.Series({
            'sample1': 1.75,
            'sample2': 3.00,
            'sample3': 1.50,
            'sample4': 2.50,
            'sample5': 4.25
        })

        with self.assertLogs('moonstone.analysis.diversity.alpha', level='WARNING') as log:
            tested_results = tested_object_instance.compute_diversity()
            pd.testing.assert_series_equal(tested_results, expected_object)
            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.analysis.diversity.alpha:INCOMPLETE TREE: missing ['species6'].\n\
Computation of the Faith's diversity using only the OTU IDs present in the Tree.", log.output)

    def test_visualize(self):
        # Not real test but make sure that visualize() runs without errors
        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 0, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tree = TreeNode.read(StringIO(
            u'(((species1:0.25,species2:0.25):0.75,species3:1.0):0.5,(species4:0.5,species5:0.5):1.0)root;'))
        tested_object_instance = FaithsPhylogeneticDiversity(tested_object, tree)
        tested_object_instance.visualize(show=False)
