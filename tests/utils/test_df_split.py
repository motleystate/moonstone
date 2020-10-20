from unittest import TestCase
import pandas as pd

from moonstone.utils.df_split import DivideByGroup


class TestDivideByGroup(TestCase):

    def test_split_into_groups(self):
        test_df = pd.DataFrame(
            [
                [23, 7, 44, 0, 101],
                [15, 4, 76, 3, 107],
                [20, 0, 22, 0, 101],
            ],
            columns=['item_1', 'item_2', 'item_3', 'item_4', 'item_5'],
            index=['1', '2', '3']  # index dtype='object'
        )

        metadata_df =  pd.DataFrame(
            [
                ['F', 1],
                ['M', 2],
                ['F', 3],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index= ['item_1', 'item_2', 'item_3', 'item_4', 'item_5'],
        )

        df1_expected = pd.DataFrame(
            [
                [23, 0],
                [15, 3],
                [20, 0]
            ],
            columns=['item_1', 'item_4'],
            index=['1', '2', '3']  # index dtype='int64'
        )

        df2_expected = pd.DataFrame(
            [
                [7, 101],
                [4, 107],
                [0, 101]
            ],
            columns=['item_2', 'item_5'],
            index=['1', '2', '3']  # index dtype='int64'
        )

        df3_expected = pd.DataFrame(
            [
                [44],
                [76],
                [22]
            ],
            columns=['item_3'],
            index=['1', '2', '3']  # index dtype='int64'
        )

        split_instance = DivideByGroup(test_df, metadata_df)
        #df1, df2, df3 = split_instance.split_df('group')
        #pd.testing.assert_frame_equal(merged_df, df_expected)