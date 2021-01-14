from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base_metaphlan import BaseMetaphlanParser


class Metaphlan2Parser(BaseMetaphlanParser):
    """
    Parse output from `Metaphlan2 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'ID'

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = self.remove_duplicates(df)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df
