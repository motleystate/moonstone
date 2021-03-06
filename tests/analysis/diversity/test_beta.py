from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.analysis.diversity.beta import (
    BrayCurtis
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
            check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
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
            [0.333, 0.714, 1.0], index=multi_index, name=BrayCurtis.DIVERSITY_INDEXES_NAME
        )
        # Two ways of retrieving the series
        pd.testing.assert_series_equal(
            tested_object_instance.beta_diversity_series, expected_object,
            check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
        )
        pd.testing.assert_series_equal(
            tested_object_instance.diversity_indexes, expected_object,
            check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
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
            check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
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
                'samples17': [8.90, 'D'],
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
        expected_object.index.names = ["Group", "Group"]

        pval = tested_object_instance._run_statistical_test_groups(
            tested_df, 'Group', stats_test='ttest_independence', correction_method='fdr_bh',
            structure_pval='series', sym=False
            )

        pd.testing.assert_series_equal(
            pval, expected_object,
            check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
        )
