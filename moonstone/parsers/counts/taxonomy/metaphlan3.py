from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base_metaphlan import BaseMetaphlanParser


class Metaphlan3Parser(BaseMetaphlanParser):
    """
    Parse output from `Metaphlan3 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'clade_name'
    NCBI_tax_column = 'NCBI_tax_id'

    def __init__(self, *args, analysis_type: str = 'rel_ab', **kwargs):
        """
        Args:
            analysis_type: output type of Metaphlan3 (see ``-t`` option of metaphlan3)
        """
        super().__init__(*args, analysis_type=analysis_type, parsing_options={'skiprows': 1}, **kwargs)

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = df.drop(self.NCBI_tax_column, axis=1)
        df = self.remove_duplicates(df)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df
