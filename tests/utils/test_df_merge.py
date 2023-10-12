from unittest import TestCase
import pandas as pd

from moonstone.utils.df_merge import MergeDF


class TestMergeDF(TestCase):

    def setUp(self):
        self.d1 = pd.DataFrame(
            [
                [23, 7, 44, 0, 101],
                [15, 4, 76, 3, 107],
                [20, 0, 22, 0, 101],
                [31, 4, 50, 0, 99]
            ],
            columns=['item_1', 'item_2', 'item_3', 'item_4', 'item_5'],
            index=[1, 2, 3, 4]  # index dtype='object'
        )

    def test_merge(self):
        d2 = pd.DataFrame(
            [
                ['M', 'Yes', 23, 'June', 170],
                ['F', 'Yes', 33, 'Nov', 154],
                ['F', 'Yes', 29, 'Jan', 161],
                ['F', 'No', 27, 'Jan', 152]
            ],
            columns=['sex', 'pets', 'age', 'sample_month', 'height'],
            index=[1, 2, 3, 4]  # index dtype='int64'
        )

        df_expected = pd.DataFrame(
            [
                ['M', 23, 7, 44, 0, 101],
                ['F', 15, 4, 76, 3, 107],
                ['F', 20, 0, 22, 0, 101],
                ['F', 31, 4, 50, 0, 99]
            ],
            columns=['sex', 'item_1', 'item_2', 'item_3', 'item_4', 'item_5'],
            index=[1, 2, 3, 4]  # index dtype='int64'
        )

        merged_df = MergeDF(self.d1, d2, 'sex').merged_df
        pd.testing.assert_frame_equal(merged_df, df_expected)

    def test_merge_index_dont_match(self):
        d2 = pd.DataFrame(
            [
                ['M', 'Yes', 23, 'June', 170],
                ['F', 'Yes', 33, 'Nov', 154],
                ['F', 'Yes', 29, 'Jan', 161],
                ['F', 'No', 27, 'Jan', 152]
            ],
            columns=['sex', 'pets', 'age', 'sample_month', 'height'],
            index=['1', '2', '3', '4']  # index dtype='object'
        )

        df_expected = pd.DataFrame(
            [
                ['M', 23, 7, 44, 0, 101],
                ['F', 15, 4, 76, 3, 107],
                ['F', 20, 0, 22, 0, 101],
                ['F', 31, 4, 50, 0, 99]
            ],
            columns=['sex', 'item_1', 'item_2', 'item_3', 'item_4', 'item_5'],
            index=['1', '2', '3', '4']  # index dtype='object'
        )

        with self.assertLogs('moonstone.utils.df_merge', level='WARNING') as log:
            merged_df = MergeDF(self.d1, d2, 'sex').merged_df
            self.assertEqual(len(log.output), 2)
            self.assertIn(
                "WARNING:moonstone.utils.df_merge:Index types do not match: int64 and object.",
                log.output
            )
        pd.testing.assert_frame_equal(merged_df, df_expected)
