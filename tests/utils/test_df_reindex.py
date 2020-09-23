from unittest import TestCase

import pandas as pd

from moonstone.utils.df_reindex import GenesToTaxonomy


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
        reindexation_instance = GenesToTaxonomy(df)  # noqa

        df_taxo = pd.DataFrame(  # noqa
            [
                [147802,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Lactobacillaceae; g__Lactobacillus; s__iners'],
                [1352,
                 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; \
f__Enterococcaceae; g__Enterococcus; s__faecium']
            ],
            columns=['full_tax', 'tax_id'],
            index=['gene_1', 'gene_2']  # index dtype='object'
        )

        # df_expected =
        # reindexed_df = reindexation_instance.reindexed_df()
        # pd.testing.assert_frame_equal(reindexed_df, df_expected)
