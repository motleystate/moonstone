from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base_metaphlan import BaseMetaphlanParser


class Metaphlan3Parser(BaseMetaphlanParser):
    """
    Parse output from `Metaphlan3 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'clade_name'
    NCBI_tax_column = 'NCBI_tax_id'

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = df.drop(self.NCBI_tax_column, axis=1)
        df = self.remove_duplicates(df)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
