from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.analysis.diversity.alpha import (
    ShannonIndex, SimpsonInverseIndex
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

        output = tested_object_instance.analyse_groups(metadata_df, 'group', make_graph=False,
                                                       output_pvalue='dataframe', sym=True)
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
                                                       output_pvalue='dataframe', sym=True,
                                                       correction_method='fdr_bh'
                                                       )
        pd.testing.assert_frame_equal(output['pval'], expected_df, check_dtype=False)


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
