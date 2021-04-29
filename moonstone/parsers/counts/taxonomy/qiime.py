from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base import BaseTaxonomyCountsParser


class Qiime2Parser(BaseTaxonomyCountsParser):
    """
    Parse output csv data obtained by `Qiime2 <https://qiime2.org/>`_.
    """

    taxa_column = '#OTU ID'
    terms_to_remove = ["Ambiguous_taxa", "Unknown Family", "uncultured"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def _load_data(self) -> DataFrame:
        df = super()._load_data()
        df = self.split_taxa_fill_none(df, sep=";", terms_to_remove=self.terms_to_remove)
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df
