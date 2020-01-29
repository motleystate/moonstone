from unittest import TestCase

import pandas as pd

from moonstone.parsers.transform.base import TransformBase


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
        col_name = 'a'
        expected_history = [
            [action, {'col_name': col_name, **arguments}]
        ]
        transform_base = TransformBase(df)
        pd.testing.assert_frame_equal(transform_base.raw_df, transform_base.df)
        self.assertFalse(transform_base.history)
        transform_base.historize(action, col_name, arguments)
        self.assertTrue(transform_base.history)
        self.assertListEqual(transform_base.history, expected_history)

    def test_rename(self):
        df = pd.DataFrame(
            [
                [1, 2, 3],
                [4, 5, 6]
            ],
            columns=['a', 'b', 'c']
        )
        expected_df = pd.DataFrame(
            [
                [1, 2, 3],
                [4, 5, 6]
            ],
            columns=['super a', 'b', 'c']
        )
        action = 'rename'
        col_name = 'a'
        arguments = {'new_name': 'super a'}
        expected_history = [
            [action, {'col_name': col_name, **arguments}]
        ]
        transform_base = TransformBase(df)
        getattr(transform_base, action)(col_name, **arguments)
        self.assertTrue(transform_base.history)
        pd.testing.assert_frame_equal(transform_base.df, expected_df)
        self.assertListEqual(transform_base.history, expected_history)
