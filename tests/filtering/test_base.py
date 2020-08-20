from unittest import TestCase

import pandas as pd

from moonstone.filtering.base import BothAxisFiltering


class FakeBothAxisFiltering(BothAxisFiltering):
    """
    Overload ABC class to test behaviour when instantiating
    """

    def filter(self):
        pass


class TestBothAxisFiltering(TestCase):

    def test_wrong_axis(self):
        test_df = pd.DataFrame.from_dict(
            {
                'specie_1': [3, 2, 1, 0],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [0, 0, 0, 0]
            },
            orient='index', columns=['1', '2', '3', '4'])
        test_df.columns.name = 'sample'
        with self.assertRaises(ValueError):
            tested_filtering = FakeBothAxisFiltering(test_df, axis=2)  # noqa
