from unittest import TestCase

import pandas as pd

from moonstone.analysis.diversity.alpha import ShannonIndex


class TestShannonIndex(TestCase):

    def test_compute_shannon_diversity(self):

        tested_object = pd.DataFrame.from_dict(
            {
                'species1': [4, 4, 0, 0, 4],
                'species2': [1, 0, 2, 0, 5],
                'species3': [0, 0, 0, 1, 4],
                'species4': [0, 3, 0, 0, 4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = ShannonIndex(tested_object)
        expected_object = pd.Series({
            'sample1': 0.721928,
            'sample2': 0.985228,
            'sample3': -0.000000,
            'sample4': -0.000000,
            'sample5': 1.992778
        })
        pd.testing.assert_series_equal(tested_object_instance.compute_alpha_diversity(), expected_object)

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
        tested_object_instance = ShannonIndex(tested_object)
        tested_object_instance.visualize(show=False)
