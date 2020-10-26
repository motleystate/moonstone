from pandas import DataFrame

from moonstone.parsers.base import BaseParser
from moonstone.utils.taxonomy import TaxonomyCountsBase


class Metaphlan2Parser(TaxonomyCountsBase, BaseParser):
    """
    Parse output from `Metaphlan2 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'ID'

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
