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
        self.analysis_type = kwargs['analysis_type']
        kwargs.pop('analysis_type')
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def remove_duplicates(self, df):
        df = df.set_index('clade_name')
        if self.analysis_type == 'rel_ab':
            total = 99.9999
        else:
            total = df[df.index.map(lambda x: len(x.split('|'))) == 1].sum()

        # dataframe at rank level
        new_df = df.index.map(lambda x: len(x.split('|')))
        self._rank_level = new_df.max()
        new_df = df[new_df == self._rank_level]

        # verification that everything is defined up to the lower_level
        samples_with_incomp_lowerlevel = new_df.sum()[new_df.sum() < total]

        rank = self._rank_level

        while samples_with_incomp_lowerlevel.size != 0:
            rank -= 1
            tmp_df_rank = df[df.index.map(lambda x: len(x.split('|'))) == rank]
            tmp_df_rank_computed = new_df.copy()
            tmp_df_rank_computed.index = tmp_df_rank_computed.index.map(lambda x: '|'.join(x.split('|')[:rank]))
            tmp_df_rank_computed = tmp_df_rank_computed.groupby(tmp_df_rank_computed.index).sum()

            rows_diff = tmp_df_rank - tmp_df_rank_computed
            rows_diff[rows_diff.isnull()] = tmp_df_rank
            rows_diff[rows_diff < 0] = 0
            rows_diff = rows_diff.loc[rows_diff.sum(axis=1)[rows_diff.sum(axis=1) != 0].index]

            if rows_diff.size != 0:
                new_df = new_df.append(rows_diff)

                # verification that everything is defined up to the lower_level
                samples_with_incomp_lowerlevel = new_df.sum()[new_df.sum() < total]

        new_df = new_df.reset_index()
        return new_df

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        df = df.drop(self.NCBI_tax_column, axis=1)
        df = self.remove_duplicates(df)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
