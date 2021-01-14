from pandas import DataFrame

from moonstone.parsers.base import BaseParser
from moonstone.utils.taxonomy import TaxonomyCountsBase


class BaseMetaphlanParser(TaxonomyCountsBase, BaseParser):

    def __init__(self, *args, analysis_type: str = 'rel_ab', **kwargs):
        """
        Args:
            analysis_type: output type of Metaphlan3 (see ``-t`` option of metaphlan3)
        """
        self.analysis_type = analysis_type
        super().__init__(*args, **kwargs)

    def rows_differences(self, dataframe1, dataframe2) -> DataFrame:
        rows_diff = dataframe1 - dataframe2
        rows_diff[rows_diff.isnull()] = dataframe1
        rows_diff[rows_diff < 0] = 0
        rows_diff = rows_diff.loc[rows_diff.sum(axis=1)[rows_diff.sum(axis=1) != 0].index]
        return rows_diff

    def compare_difference_between_two_levels(self, whole_df, df_at_lower_level, rank) -> DataFrame:
        df_rank = whole_df[whole_df.index.map(lambda x: len(x.split('|'))) == rank]

        # transformation lower_level to rank (level)
        df_rank_computed = df_at_lower_level.copy()
        df_rank_computed.index = df_rank_computed.index.map(lambda x: '|'.join(x.split('|')[:rank]))   # to rank (level)
        df_rank_computed = df_rank_computed.groupby(df_rank_computed.index).sum()             # grouping by rank (level)
        return self.rows_differences(df_rank, df_rank_computed)

    def remove_duplicates(self, df) -> DataFrame:
        df = df.set_index(self.taxa_column)

        # dataframe at rank level
        index_levels = df.index.map(lambda x: len(x.split('|')))          # first, creation of the index
        self.rank_level = index_levels.max()                              # max rank level
        first_rank = index_levels.min()
        new_df = df[index_levels == self.rank_level]

        # calculation of the total
        if self.analysis_type == 'rel_ab':
            total = 99.9999                     # addition error margin
        else:
            total = df[df.index.map(lambda x: len(x.split('|'))) == first_rank].sum()

        # verification that everything is defined up to the lower_level
        samples_with_incomp_lowerlevel = new_df.sum()[new_df.sum() < total]

        rank = self.rank_level

        while samples_with_incomp_lowerlevel.size != 0 and rank > 1:
            rank -= 1
            rows_diff = self.compare_difference_between_two_levels(df, new_df, rank)
            if rows_diff.size != 0:
                new_df = new_df.append(rows_diff)              # add missing rows to the dataframe of the lower level

                # verification that everything is defined up to the lower_level
                samples_with_incomp_lowerlevel = new_df.sum()[new_df.sum() < total]

        if samples_with_incomp_lowerlevel.size != 0:      # for rel ab, if < 100 at the kingdom level
            rows_diff = 100 - new_df.sum()
            rows_diff.name = 'k__UNKNOWN'
            new_df = new_df.append(rows_diff)

        new_df = new_df.reset_index()
        return new_df
