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

    def test_invalid_pval_param(self):
        tested_object_instance = ShannonIndex(self.tested_object)
        with self.assertLogs('moonstone.analysis.diversity.base', level='WARNING') as log:
            tested_object_instance._valid_pval_param("lalala")
            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.analysis.diversity.base:pval_to_compute='lalala' not valid, \
set to default (all).", log.output)

    def test_analyse_groups_pval_to_compute_all(self):
        tested_object_instance = ShannonIndex(self.tested_object)

        if version.parse(scipy.__version__) > version.parse("1.5.3"):
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
        # if scipy 1.5.3 = no `method` argument -> uses "asymptotic" method
        # if scipy > 1.5.3 (ex: 1.7.1) = `method` argument; default being 'auto',
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
        if version.parse(scipy.__version__) > version.parse("1.5.3"):
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
        # if scipy 1.5.3 = no `method` argument -> uses "asymptotic" method
        # if scipy > 1.5.3 (ex: 1.7.1) = `method` argument; default being 'auto',
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
