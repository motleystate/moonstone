import pandas as pd

from moonstone.parsers.counts.taxonomy.base import TaxonomyCountsBaseParser


class SunbeamKraken2Parser(TaxonomyCountsBaseParser):
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
