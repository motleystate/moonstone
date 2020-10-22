from unittest import TestCase

import pandas as pd

from moonstone.analysis.diversity.beta import (
    BrayCurtis
)


class TestBrayCurtis(TestCase):

    def test_compute_beta_diversity_df(self):

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0],
                'species2': [1, 0, 2],
                'species3': [0, 0, 0],
                'species4': [0, 3, 0]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3'])
        tested_object_instance = BrayCurtis(tested_object)
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

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0],
                'species2': [1, 0, 2],
                'species3': [0, 0, 0],
                'species4': [0, 3, 0]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3'])
        tested_object_instance = BrayCurtis(tested_object)
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
        pd.testing.assert_series_equal(
            tested_object_instance.beta_diversity_series, expected_object,
            check_less_precise=2,  # Deprecated since version 1.1.0, to be changed when updating pandas
        )
