from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.utils.pandas.df_reindex import GenesToTaxonomy


class TestGenesToTaxonomy(TestCase):

    def test_reindex_with_taxonomy(self):
        df = pd.DataFrame(
            [
                [23, 7],
                [15, 4],
            ],
            columns=['sample_1', 'sample_2'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_taxo = pd.DataFrame(
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [1352,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Enterococcaceae; g__Enterococcus; s__faecium']
            ],
            columns=['tax_id', 'full_tax'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_expected = pd.DataFrame.from_dict(
            {
                'sample_1':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 15,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 23
                },
                'sample_2':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Enterococcaceae', 'Enterococcus', 'Enterococcus_faecium'): 4,
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 7
                }
            }
        )
        df_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        reindexation_instance = GenesToTaxonomy(df, df_taxo)
        reindexed_df = reindexation_instance.reindexed_df
        pd.testing.assert_frame_equal(reindexed_df, df_expected)

    def test_reindex_with_taxonomy_missing_infos_dropped(self):
        df = pd.DataFrame(
            [
                [23, 7],
                [15, 4],
            ],
            columns=['sample_1', 'sample_2'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_taxo = pd.DataFrame(
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [1352,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Enterococcaceae; g__Enterococcus; s__faecium']
            ],
            columns=['tax_id', 'full_tax'],
            index=['gene_1', 'gene_3']  # index dtype='object'
        )
        df_expected = pd.DataFrame.from_dict(
            {
                'sample_1':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 23
                },
                'sample_2':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 7
                }
            }
        )
        df_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        reindexation_instance = GenesToTaxonomy(df, df_taxo)
        reindexed_df = reindexation_instance.reindexed_df
        pd.testing.assert_frame_equal(reindexed_df, df_expected)
        pd.testing.assert_index_equal(reindexation_instance.without_info_index, pd.Index(['gene_2'], dtype='object'))

    def test_reindex_with_taxonomy_missing_infos_kept(self):
        df = pd.DataFrame(
            [
                [23, 7],
                [15, 4],
                [0, 36],
            ],
            columns=['sample_1', 'sample_2'],
            index=['gene_1', 'gene_2', 'gene_4']  # index dtype='object'
        )
        df_taxo = pd.DataFrame(
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [1352,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Enterococcaceae; g__Enterococcus; s__faecium']
            ],
            columns=['tax_id', 'full_tax'],
            index=['gene_1', 'gene_3']  # index dtype='object'
        )
        df_expected = pd.DataFrame.from_dict(
            {
                'sample_1':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 23,
                    (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'gene_2_species'): 15,
                    (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'gene_4_species'): 0,
                },
                'sample_2':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 7,
                    (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'gene_2_species'): 4,
                    (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'gene_4_species'): 36,
                }
            }
        )
        df_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        reindexation_instance = GenesToTaxonomy(df, df_taxo)
        reindexed_df = reindexation_instance.reindex_with_taxonomy(na='keep')
        pd.testing.assert_frame_equal(reindexed_df, df_expected)
        pd.testing.assert_index_equal(
            reindexation_instance.without_info_index,
            pd.Index(['gene_2', 'gene_4'], dtype='object')
        )

    def test_reindex_with_taxonomy_missing_infos_summed(self):
        df = pd.DataFrame(
            [
                [23, 7],
                [15, 4],
                [0, 36],
            ],
            columns=['sample_1', 'sample_2'],
            index=['gene_1', 'gene_2', 'gene_4']  # index dtype='object'
        )
        df_taxo = pd.DataFrame(
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [1352,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Enterococcaceae; g__Enterococcus; s__faecium']
            ],
            columns=['tax_id', 'full_tax'],
            index=['gene_1', 'gene_3']  # index dtype='object'
        )
        df_expected = pd.DataFrame.from_dict(
            {
                'sample_1':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 23,
                    (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan): 15,
                },
                'sample_2':
                {
                    ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                     'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 7,
                    (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan): 40,
                }
            }
        )
        df_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        reindexation_instance = GenesToTaxonomy(df, df_taxo)
        reindexed_df = reindexation_instance.reindex_with_taxonomy(na='sum')
        pd.testing.assert_frame_equal(reindexed_df, df_expected)

    def test_reindex_with_taxonomy_summing(self):
        df = pd.DataFrame(
            [
                [23, 7],
                [15, 4],
            ],
            columns=['sample_1', 'sample_2'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_taxo = pd.DataFrame(
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners']
            ],
            columns=['tax_id', 'full_tax'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_expected = pd.DataFrame.from_dict(
             {
                 'sample_1':
                 {
                     ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                      'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 38
                 },
                 'sample_2':
                 {
                     ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                      'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 11,
                 }
             }
         )
        df_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        reindexation_instance = GenesToTaxonomy(df, df_taxo)
        reindexed_df = reindexation_instance.reindexed_df
        pd.testing.assert_frame_equal(reindexed_df, df_expected)

    def test_reindex_with_taxonomy_counting(self):
        df = pd.DataFrame(
            [
                [23, 7],
                [15, 0],
            ],
            columns=['sample_1', 'sample_2'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_taxo = pd.DataFrame(
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners']
            ],
            columns=['tax_id', 'full_tax'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )
        df_expected = pd.DataFrame.from_dict(
             {
                 'sample_1':
                 {
                     ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                      'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 2
                 },
                 'sample_2':
                 {
                     ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales',
                      'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners'): 1,
                 }
             }
         )
        df_expected.index.set_names(["kingdom", "phylum", "class", "order", "family", "genus", "species"], inplace=True)

        reindexation_instance = GenesToTaxonomy(df, df_taxo)
        reindexed_df = reindexation_instance.reindex_with_taxonomy(method='count')
        pd.testing.assert_frame_equal(reindexed_df, df_expected)
