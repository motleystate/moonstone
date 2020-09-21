import pandas as pd

from moonstone.utils.taxonomy import TaxonomyCountsBase


class GenesToTaxonomy(TaxonomyCountsBase):

    def __init__(self, dataframe):
        self.df = dataframe
        taxa_column = 'full_tax'

    def reindex_with_taxonomy(self, taxonomy_file):
        df_taxonomy = pd.read_csv(taxonomy_file, index_col=0)
        new_df = self.df.merge(df_taxonomy['full_tax'], how='left', left_index=True, right_index=True)
        new_df['full_tax'] = new_df['full_tax'].fillna(value='k__; p__; c__; o__; f__; g__; s__')
        new_df = self.split_taxa_fill_none(new_df, sep="; ", merge_genus_species=True)
        new_df = new_df.set_index(self.taxonomical_names[:self._rank_level])