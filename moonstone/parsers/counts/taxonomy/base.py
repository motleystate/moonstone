from typing import List

import numpy as np
import pandas as pd

from moonstone.parsers.base import BaseParser


class TaxonomyCountsBaseParser(BaseParser):

    taxonomical_names = [
        "kingdom", "phylum", "class", "order", "family", "genus", "species", "sTrain"
    ]
    taxa_column = 'OTU ID'

    def _fill_none(self, taxa_df: pd.DataFrame) -> pd.DataFrame:
        """
        This method is used to obtain a dataframe that fills the None values with the last valid value and
        the category where it belonged to.

        e.g:

        Before:
        column names     kingdom  phylum            family         genus
        value            Bacteria Bacteroidetes ... Tannerellaceae None

        After:
        column names     kingdom  phylum            family         genus
        value            Bacteria Bacteroidetes ... Tannerellaceae Tannerellaceae (family)
        """
        taxa_df_with_rank = taxa_df.apply(lambda x: x + " ({})".format(x.name))
        taxa_df_with_rank_filled_none = taxa_df_with_rank.fillna(method='ffill', axis=1)
        taxa_df_filled_none = taxa_df.combine_first(taxa_df_with_rank_filled_none)
        return taxa_df_filled_none

    def _merge_genus_species(self, taxa_df: pd.DataFrame) -> pd.DataFrame:

        def join_genus_species(genus_and_species: List):
            """
            :param genus_and_species: list of 2 items
            """
            if genus_and_species[1] is None:
                return None
            elif genus_and_species[0] is None:
                return genus_and_species[1]
            else:
                return '_'.join(genus_and_species)

        taxa_df['species'] = taxa_df[['genus', 'species']].apply(lambda x: join_genus_species(x), axis=1)
        return taxa_df

    def split_taxa_fill_none(self, df: pd.DataFrame, sep: str = ";",  taxo_prefix: str = "__",
                             merge_genus_species: bool = False, terms_to_remove: List = None) -> pd.DataFrame:
        """
        :param terms_to_remove: if specified, list of term to remove from taxa names (e.g. uncultured)
        """

        def remove_taxo_prefix(string):
            if string is None:
                return None
            else:
                term = string.split(taxo_prefix)[-1]
            if term:
                return term
            else:
                return None

        taxa_columns = df[self.taxa_column].str.split(sep, expand=True)
        self._rank_level = len(taxa_columns.columns)
        taxa_columns.columns = self.taxonomical_names[:self._rank_level]
        taxa_columns = taxa_columns.applymap(lambda x: remove_taxo_prefix(x))
        if terms_to_remove is not None:
            taxa_columns = taxa_columns.replace(terms_to_remove, np.nan)
        if merge_genus_species:
            taxa_columns = self._merge_genus_species(taxa_columns)
        taxa_columns = self._fill_none(taxa_columns)
        return pd.concat([self._fill_none(taxa_columns), df.drop(self.taxa_column, axis=1)], axis=1)
