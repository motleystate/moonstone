from pandas import DataFrame

from moonstone.parsers.base import BaseParser
from moonstone.utils.taxonomy import TaxonomyCountsBase


class Metaphlan3Parser(TaxonomyCountsBase, BaseParser):
    """
    Parse output from `Metaphlan3 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'clade_name'
    NCBI_tax_column = 'NCBI_tax_id'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = df.drop(self.NCBI_tax_column, axis=1)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
