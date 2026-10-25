import os
from unittest import TestCase

import pandas as pd

from moonstone.parsers.counts.taxonomy import MetaphlanReadStatsParser

DIR = os.path.dirname(__file__)
LEVELS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'sTrain']
LACTO = ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus']
GARD = ['Bacteria', 'Actinobacteria', 'Actinomycetia', 'Bifidobacteriales', 'Bifidobacteriaceae', 'Gardnerella']
UNCLASSIFIED = ['UNCLASSIFIED'] + ['UNCLASSIFIED (kingdom)'] * 7


class TestMetaphlanReadStatsParser(TestCase):

    def setUp(self):
        self.path_a = os.path.join(DIR, 'sample_A_profile.txt')
        self.path_b = os.path.join(DIR, 'sample_B_profile.txt')

    def test_to_dataframe(self):
        """
        Two MetaPhlAn 4 profiles: species are joined on clade name, a clade absent from a sample has 0 reads,
        and the UNCLASSIFIED reads are kept, so each sample sums to its total number of reads.
        """
        parser = MetaphlanReadStatsParser([self.path_a, self.path_b])
        expected_df = pd.DataFrame(
            [
                LACTO + ['Lactobacillus_crispatus', 'SGB7037', 6000, 0],
                LACTO + ['Lactobacillus_iners', 'SGB7046', 2000, 1000],
                GARD + ['Gardnerella_vaginalis', 'SGB17302', 1000, 4000],
                UNCLASSIFIED + [1000, 0],
            ],
            columns=LEVELS + ['sample_A', 'sample_B']
        ).set_index(LEVELS)
        pd.testing.assert_frame_equal(parser.dataframe, expected_df, check_like=True)
        self.assertEqual(parser.dataframe.sum().to_dict(), {'sample_A': 10000, 'sample_B': 5000})

    def test_sample_names(self):
        parser = MetaphlanReadStatsParser([self.path_a, self.path_b], sample_names=['IMVM_1', 'IMVM_2'])
        self.assertEqual(list(parser.dataframe.columns), ['IMVM_1', 'IMVM_2'])

    def test_single_file_path(self):
        parser = MetaphlanReadStatsParser(self.path_b)
        self.assertEqual(list(parser.dataframe.columns), ['sample_B'])
        self.assertEqual(parser.dataframe['sample_B'].sum(), 5000)

    def test_estimated_bases(self):
        """
        With long reads, MetaPhlAn estimates bases instead of reads.
        """
        parser = MetaphlanReadStatsParser(os.path.join(DIR, 'sample_C_long_reads.txt'))
        expected_df = pd.DataFrame(
            [
                LACTO + ['Lactobacillus_crispatus', 'SGB7037', 600],
                UNCLASSIFIED + [200],
            ],
            columns=LEVELS + ['sample_C_long_reads']
        ).set_index(LEVELS)
        pd.testing.assert_frame_equal(parser.dataframe, expected_df, check_like=True)

    def test_fewer_taxonomical_names(self):
        """
        Rows deeper than taxonomical_names (here the SGB level) are left out.
        """
        class SpeciesParser(MetaphlanReadStatsParser):
            taxonomical_names = LEVELS[:7]

        observed_df = SpeciesParser([self.path_a, self.path_b]).dataframe
        self.assertEqual(list(observed_df.index.names), LEVELS[:7])
        self.assertEqual(observed_df.loc[tuple(LACTO + ['Lactobacillus_iners'])].tolist(), [2000, 1000])
        self.assertEqual(observed_df.sum().to_dict(), {'sample_A': 10000, 'sample_B': 5000})

    def test_wrong_number_of_sample_names(self):
        with self.assertRaisesRegex(ValueError, "1 sample_names given for 2 files"):
            MetaphlanReadStatsParser([self.path_a, self.path_b], sample_names=['IMVM_1'])

    def test_duplicated_sample_names(self):
        with self.assertRaisesRegex(ValueError, "Duplicated sample names"):
            MetaphlanReadStatsParser([self.path_a, self.path_b], sample_names=['IMVM_1', 'IMVM_1'])

    def test_not_a_read_stats_profile(self):
        parser = MetaphlanReadStatsParser(os.path.join(DIR, 'rel_ab_profile.txt'))
        with self.assertRaisesRegex(ValueError, "is not a MetaPhlAn rel_ab_w_read_stats profile"):
            parser.dataframe

    def test_empty_profile(self):
        parser = MetaphlanReadStatsParser(os.path.join(DIR, 'empty_profile.txt'))
        with self.assertRaisesRegex(ValueError, "has no clade rows"):
            parser.dataframe
