from unittest import TestCase

import pandas as pd
import numpy as np

from moonstone.filtering.taxonomy_filtering import TaxonomyNamesFiltering


class TestTaxonomyNamesFiltering(TestCase):

    def test_selecting_rows_multiindex_species_level(self):
        test_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        test_df.columns.name = 'sample'
        test_df.index.set_names(['genus', 'species'], inplace=True)
        selected_rows = ['specie_1', 'specie_2']
        tested_filtering = TaxonomyNamesFiltering(test_df, selected_rows, level='species', keep=True)
        expected_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
            },
            orient='index')
        expected_df.columns.name = 'sample'
        expected_df.index.set_names(['genus', 'species'], inplace=True)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_exclude_rows_multiindex_species_level(self):
        test_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        test_df.columns.name = 'sample'
        test_df.index.set_names(['genus', 'species'], inplace=True)
        selected_rows = ['specie_1', 'specie_2']
        tested_filtering = TaxonomyNamesFiltering(test_df, selected_rows, level='species', keep=False)
        expected_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        expected_df.columns.name = 'sample'
        expected_df.index.set_names(['genus', 'species'], inplace=True)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_selecting_rows_multiindex_genus_level(self):
        test_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        test_df.columns.name = 'sample'
        test_df.index.set_names(['genus', 'species'], inplace=True)
        selected_rows = ['genus_1']
        tested_filtering = TaxonomyNamesFiltering(test_df, selected_rows, level='genus', keep=True)
        expected_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
            },
            orient='index')
        expected_df.columns.name = 'sample'
        expected_df.index.set_names(['genus', 'species'], inplace=True)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_exclude_rows_multiindex_genus_level(self):
        test_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        test_df.columns.name = 'sample'
        test_df.index.set_names(['genus', 'species'], inplace=True)
        selected_rows = ['genus_1']
        tested_filtering = TaxonomyNamesFiltering(test_df, selected_rows, level='genus', keep=False)
        expected_df = pd.DataFrame.from_dict(
            {
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        expected_df.columns.name = 'sample'
        expected_df.index.set_names(['genus', 'species'], inplace=True)
        pd.testing.assert_frame_equal(tested_filtering.filtered_df, expected_df)

    def test_wrong_level(self):
        test_df = pd.DataFrame.from_dict(
            {
                ('genus_1', 'specie_1'): {'1': 3, '2': 2, '3': 1, '4': 0},
                ('genus_1', 'specie_3'): {'1': 5, '2': 8, '3': 4, '4': 7},
                ('genus_2', 'specie_2'): {'1': 25, '2': 6, '3': 3, '4': 9},
                ('genus_2', 'specie_4'): {'1': 9, '2': 18, '3': 6, '4': 3},
                (np.nan, np.nan): {'1': 4, '2': 8, '3': 9, '4': 0}
            },
            orient='index')
        test_df.columns.name = 'sample'
        test_df.index.set_names(['genus', 'species'], inplace=True)
        selected_rows = ['genus_1']
        with self.assertRaises(ValueError):
            tested_filtering = TaxonomyNamesFiltering(test_df, selected_rows, level='phylum', keep=False)  # noqa
