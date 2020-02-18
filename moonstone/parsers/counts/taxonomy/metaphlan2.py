import pandas as pd

from moonstone.parsers.counts.taxonomy.base import TaxonomyCountsBaseParser


class Metaphlan2Parser(TaxonomyCountsBaseParser):
    """
    Parse output from metaphlan2 merged table
    """

    taxonomical_names = [
        "kingdom", "phylum", "class", "order", "family", "genus", "species", "sTrain"
    ]
    taxa_column = 'ID'

    def split_taxa_fill_none(self, df):
        """
        This function split taxa column into different ones.
        It also fill in None with latest found annotation
        """
        def remove_taxo_prefix(string):
            if string is None:
                return None
            else:
                return string.split('__')[-1]

        taxa_columns = df[self.taxa_column].str.split("|", expand=True)
        self._rank_level = len(taxa_columns.columns)
        taxa_columns.columns = self.taxonomical_names[:self._rank_level]
        taxa_columns = taxa_columns.applymap(lambda x: remove_taxo_prefix(x))
        taxa_columns = self._fill_none(taxa_columns)
        return pd.concat([self._fill_none(taxa_columns), df.drop(self.taxa_column, axis=1)], axis=1)

    def to_dataframe(self):
        df = super().to_dataframe()
        df = self.split_taxa_fill_none(df)
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
