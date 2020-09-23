import pandas as pd
from typing import Union

from moonstone.utils.taxonomy import TaxonomyCountsBase


class GenesToTaxonomy(TaxonomyCountsBase):

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame], taxonomy_file: str = None):
        self.df = dataframe
        if taxonomy_file is not None:
            self.taxonomy_file = taxonomy_file
        taxa_column = 'full_tax'   # noqa

    def reindex_with_taxonomy(self):
        new_df = self.df.merge(self.taxonomy_df['full_tax'], how='left', left_index=True, right_index=True)
        new_df['full_tax'] = new_df['full_tax'].fillna(value='k__; p__; c__; o__; f__; g__; s__')
        new_df = self.split_taxa_fill_none(new_df, sep="; ", merge_genus_species=True)
        new_df = new_df.set_index(self.taxonomical_names[:self._rank_level])
        return new_df

    @property
    def reindexed_df(self):
        if getattr(self, "_reindexed_df", None) is None:
            self._reindexed_df = self.reindex_with_taxonomy()
        return self._reindexed_df

    @property
    def taxonomy_df(self):
        """
        retrieves taxonomy_df, and read it from given taxonomy_file if no values given
        """
        if getattr(self, "_taxonomy_df", None) is None:
            self._taxonomy_df = pd.read_csv(self.taxonomy_file, index_col=0)
        return self._taxonomy_df

    @taxonomy_df.setter
    def taxonomy_df(self, taxonomy_df):
        # if blabla:
        #    # check that index correspond to df's index
        # else:
        #    raise ValueError("Error : given taxonomy_df's index doesn't match df's index")
        pass
