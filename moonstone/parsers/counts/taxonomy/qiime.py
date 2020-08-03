from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base import TaxonomyCountsBaseParser


class Qiime2Parser(TaxonomyCountsBaseParser):
    """
    This class is used to parse the data obtained by Qiime2 (https://qiime2.org/).
    """

    taxa_column = '#OTU ID'
    terms_to_remove = ["Ambiguous_taxa", "Unknown Family", "uncultured"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = self.split_taxa_fill_none(df, sep=";", terms_to_remove=self.terms_to_remove)
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
