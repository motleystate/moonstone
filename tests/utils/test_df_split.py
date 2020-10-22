from unittest import TestCase
import pandas as pd

from moonstone.utils.df_split import DivideByGroup


class TestDivideByGroup(TestCase):

    def setUp(self):
        self.metadata_df = pd.DataFrame(
            [
                ['F', 1],
                ['M', 2],
                ['F', 3],
                ['M', 1],
                ['M', 2]
            ],
            columns=['sex', 'group'],
            index=['sample_1', 'sample_2', 'sample_3', 'sample_4', 'sample_5'],
        )

    def test_one_group_one_dataframe(self):
        test_df = pd.DataFrame(
            [
                [23, 7, 44, 0, 101],
                [15, 4, 76, 3, 107],
                [20, 0, 22, 0, 101],
            ],
            columns=['sample_1', 'sample_2', 'sample_3', 'sample_4', 'sample_5'],
            index=['1', '2', '3']  # index dtype='object'
        )

        df1_expected = pd.DataFrame(
            [
                [23, 0],
                [15, 3],
                [20, 0]
            ],
            columns=['sample_1', 'sample_4'],
            index=['1', '2', '3']  # index dtype='int64'
        )
        df2_expected = pd.DataFrame(
            [
                [7, 101],
                [4, 107],
                [0, 101]
            ],
            columns=['sample_2', 'sample_5'],
            index=['1', '2', '3']  # index dtype='int64'
        )
        df3_expected = pd.DataFrame(
            [
                [44],
                [76],
                [22]
            ],
            columns=['sample_3'],
            index=['1', '2', '3']  # index dtype='int64'
        )

        split_instance = DivideByGroup(test_df, self.metadata_df)
        df1, df2, df3 = split_instance.split_df('group')
        pd.testing.assert_frame_equal(df1, df1_expected)
        pd.testing.assert_frame_equal(df2, df2_expected)
        pd.testing.assert_frame_equal(df3, df3_expected)

    def test_more_than_one_group_in_dataframe(self):
        test_df = pd.DataFrame(
            [
                [23, 7, 44, 0, 101],
                [15, 4, 76, 3, 107],
                [20, 0, 22, 0, 101],
            ],
            columns=['sample_1', 'sample_2', 'sample_3', 'sample_4', 'sample_5'],
            index=['1', '2', '3']  # index dtype='object'
        )

        df1_expected = pd.DataFrame(
            [
                [23, 0],
                [15, 3],
                [20, 0]
            ],
            columns=['sample_1', 'sample_4'],
            index=['1', '2', '3']  # index dtype='int64'
        )
        df23_expected = pd.DataFrame(
            [
                [7, 44, 101],
                [4, 76, 107],
                [0, 22, 101]
            ],
            columns=['sample_2', 'sample_3', 'sample_5'],
            index=['1', '2', '3']  # index dtype='int64'
        )

        split_instance = DivideByGroup(test_df, self.metadata_df)
        df1, df23 = split_instance.split_df('group', division_seq='1_2-3')
        pd.testing.assert_frame_equal(df1, df1_expected)
        pd.testing.assert_frame_equal(df23, df23_expected)

    def test_taxonomy_dataframe(self):
        test_df = pd.DataFrame.from_dict(
            {
                'sample_1':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 23,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 15
                },
                'sample_2':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 7,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 4
                },
                'sample_3':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 98,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 49
                },
                'sample_4':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 17,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 0
                },
                'sample_5':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 21,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 3
                }
            }
        )
        test_df.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        df1_expected = pd.DataFrame.from_dict(
            {
                'sample_1':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 23,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 15
                },
                'sample_4':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 17,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 0
                }
            }
        )
        df1_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"],
                                     inplace=True)
        df2_expected = pd.DataFrame.from_dict(
            {
                'sample_2':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 7,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 4
                },
                'sample_5':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 21,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 3
                }
            }
        )
        df2_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"],
                                     inplace=True)
        df3_expected = pd.DataFrame.from_dict(
            {
                'sample_3':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 98,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 49
                }
            }
        )
        df3_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"],
                                     inplace=True)

        split_instance = DivideByGroup(test_df, self.metadata_df)
        df1, df2, df3 = split_instance.split_df('group')
        pd.testing.assert_frame_equal(df1, df1_expected, check_dtype=False)
        pd.testing.assert_frame_equal(df2, df2_expected, check_dtype=False)
        pd.testing.assert_frame_equal(df3, df3_expected, check_dtype=False, check_like=True)  # ignore order of index
