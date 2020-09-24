import pandas as pd
from typing import Union

from moonstone.utils.taxonomy import TaxonomyCountsBase


class GenesToTaxonomy(TaxonomyCountsBase):

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame],
                 taxonomy_file: str = None, taxa_column: str = 'full_tax'):
        self.df = dataframe
        self.taxonomy_file = taxonomy_file
        self.taxa_column = taxa_column   # noqa

    def stats_on_taxonomy_merge(self, merged_df):
        """
        :param merged_df needs to have been merged with the indicator arg set to True
        """
        tot = merged_df['_merge'].size
        both_counts = merged_df['_merge'].value_counts()['both']
        print(f"Merge of taxonomic data to the count dataframe for {both_counts} items ({(both_counts/tot)*100}%).")
        if both_counts != tot:
            print("If these results aren't as expected, \
please check that the indexes (items name) of both dataframes match.")
            print("You can access the list of items without taxonomic infos \
by checking the .without_infos_index attribute.")
            self.without_infos_index = merged_df.loc[merged_df['_merge'] == 'left_only'].index

    def reindex_with_taxonomy(self, stats: bool = True):
        """
        :param stats: boolean to show merging stats
        """
        new_df = self.df.merge(self.taxonomy_df[self.taxa_column], how='left',
                               left_index=True, right_index=True, indicator=stats)
        if stats is True:
            self.stats_on_taxonomy_merge(new_df)
            new_df = new_df.drop(['_merge'], axis=1)
        new_df[self.taxa_column] = new_df[self.taxa_column].fillna(value='k__; p__; c__; o__; f__; g__; s__')
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
        retrieves taxonomy_df, and read it from given taxonomy_file if no value has been given
        """
        if getattr(self, "_taxonomy_df", None) is None:
            if self.taxonomy_file is not None:
                self._taxonomy_df = pd.read_csv(self.taxonomy_file)
                self._taxonomy_df = self._taxonomy_df.set_index(self.taxa_column)
            else:
                raise ValueError("No taxonomy_df nor taxonomy_file has been specified")
        return self._taxonomy_df

    @taxonomy_df.setter
    def taxonomy_df(self, taxonomy_df):
        self._taxonomy_df = taxonomy_df
