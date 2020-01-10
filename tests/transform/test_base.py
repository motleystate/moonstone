from unittest import TestCase

import pandas as pd

from moonstone.transform.base import TransformBase


class TestTransformBase(TestCase):

    def test_historize(self):
        df = pd.DataFrame(
            [
                [1, 2, 3],
                [4, 5, 6]
            ],
            columns=['a', 'b', 'c']
        )
        action = 'an_action'
        arguments = {'arg1': 1, 'arg2': 2}
        expected_history = [
            [action, arguments]
        ]
        transform_base = TransformBase(df)
        pd.testing.assert_frame_equal(transform_base.raw_df, transform_base.df)
        self.assertFalse(transform_base.history)
        transform_base.historize(action, arguments)
        self.assertTrue(transform_base.history)
        self.assertListEqual(transform_base.history, expected_history)
