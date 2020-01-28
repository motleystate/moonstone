from unittest import TestCase

import pandas as pd

from moonstone.parsers.transform.cleaning import StringCleaner


class TestStringCleaner(TestCase):

    def test_remove_trailing_spaces(self):
        df = pd.DataFrame(
            [
                [1, ' b'],
                [4, " a "]
            ],
            columns=['number', 'string']
        )
        expected_df = pd.DataFrame(
            [
                [1, 'b'],
                [4, "a"]
            ],
            columns=['number', 'string']
        )
        method_name = "remove_trailing_spaces"
        expected_history = [
            [method_name, {'col_name': 'string'}]
        ]
        transform_cleaning = StringCleaner(df)
        getattr(transform_cleaning, method_name)('string')
        self.assertTrue(transform_cleaning.history)
        self.assertListEqual(transform_cleaning.history, expected_history)
        pd.testing.assert_frame_equal(transform_cleaning.df, expected_df)

    def test_to_slug(self):
        df = pd.DataFrame(
            [
                [1, ' b test  '],
                [4, " a Stuff.2"]
            ],
            columns=['number', 'string']
        )
        expected_df = pd.DataFrame(
            [
                [1, 'b-test'],
                [4, "a-stuff-2"]
            ],
            columns=['number', 'string']
        )
        method_name = "to_slug"
        expected_history = [
            [method_name, {'col_name': 'string'}]
        ]
        transform_cleaning = StringCleaner(df)
        getattr(transform_cleaning, method_name)('string')
        self.assertTrue(transform_cleaning.history)
        self.assertListEqual(transform_cleaning.history, expected_history)
        pd.testing.assert_frame_equal(transform_cleaning.df, expected_df)
