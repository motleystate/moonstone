from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base import TaxonomyCountsBaseParser


class Metaphlan2Parser(TaxonomyCountsBaseParser):
    """
    Parse output from metaphlan2 merged table
    """

    taxa_column = 'ID'

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
