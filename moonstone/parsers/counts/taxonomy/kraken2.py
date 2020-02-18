import pandas as pd
from moonstone.parsers.base import BaseParser


class SunbeamKraken2Parser(BaseParser):
    """
    Parse output from Kraken2 merge table from Sunbeam pipeline
    """

    taxonomical_names = [
        "kingdom", "phylum", "class", "order", "family", "genus", "species"
    ]
    taxa_column = 'Consensus Lineage'
    new_otu_id_name = 'NCBI_taxonomy_ID'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def _fill_none(self, taxa_df):
        """
        This function serves to obtain a data frame that fills the None values with the last valid value and
        the category where it belonged to. E.g:

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

    def split_taxa_fill_none(self, df):
        """
        This function split taxa column into different ones.
        It also fill in None with latest found annotation
        """
        def remove_taxo_prefix(string):
            if string is None:
                return None
            else:
                term = string.split('__')[-1]
            if term:
                return term
            else:
                return None

        taxa_columns = df[self.taxa_column].str.split("; ", expand=True)
        self._rank_level = len(taxa_columns.columns)
        taxa_columns.columns = self.taxonomical_names[:self._rank_level]
        taxa_columns = taxa_columns.applymap(lambda x: remove_taxo_prefix(x))
        taxa_columns = self._fill_none(taxa_columns)
        return pd.concat([self._fill_none(taxa_columns), df.drop(self.taxa_column, axis=1)], axis=1)

    def to_dataframe(self):
        df = super().to_dataframe()
        # Rename first column to NCBI_taxonomy_ID
        df.columns = [self.new_otu_id_name] + list(df.columns[1:])
        df = self.split_taxa_fill_none(df)
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
