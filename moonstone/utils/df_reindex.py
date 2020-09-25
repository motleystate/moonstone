import logging
import pandas as pd
from typing import Union

from moonstone.utils.taxonomy import TaxonomyCountsBase

logger = logging.getLogger(__name__)



class GenesToTaxonomy(TaxonomyCountsBase):

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame],
                 taxonomy_dataframe: Union[pd.Series, pd.DataFrame], taxa_column: str = 'full_tax'):
        """
        :param dataframe : genes counts dataframe
        :param taxonomy_dataframe : index items (genes) names and with a column with the full taxonomic informations
        following the `Kraken2 <https://ccb.jhu.edu/software/kraken2/>`_ format 'k__; p__; c__; o__; f__; g__; s__
        :param taxa_column : name of the column in taxonomy_dataframe containing the full taxonomic informations
        """
        self.df = dataframe
        self.taxonomy_df = taxonomy_dataframe
        self.taxa_column = taxa_column

    def reindex_with_taxonomy(self):
        """
        reindexation on taxonomic informations (if there are).

        NB: You can access the list of items without taxonomic infos by checking the .without_infos_index attributes
        """
        # merging of count dataframe with taxonomy dataframe
        new_df = self.df.merge(self.taxonomy_df[self.taxa_column], how='left',
                               left_index=True, right_index=True, indicator=True)
        # stats and warnings on merge
        tot = new_df['_merge'].size
        both_counts = new_df['_merge'].value_counts()['both']
        logger.info(f"Merge of taxonomic data to the count dataframe for {both_counts} items ({(both_counts/tot)*100}%).")
        if both_counts != tot:
            logger.info("If these results aren't as expected, \
please check that the indexes (items name) of both dataframes match.")
            logger.info("You can access the list of items without taxonomic infos \
by checking the .without_infos_index attribute.")
        self.without_infos_index = new_df['_merge'].loc[new_df['_merge'] == 'left_only'].index

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
