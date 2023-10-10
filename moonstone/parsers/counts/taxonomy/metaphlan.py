import logging

import pandas as pd

from moonstone.parsers.counts.taxonomy.base import BaseTaxonomyCountsParser

logger = logging.getLogger(__name__)


class BaseMetaphlanParser(BaseTaxonomyCountsParser):

    def __init__(self, *args, analysis_type: str = 'rel_ab', **kwargs):
        """
        Args:
            analysis_type: output type of Metaphlan3 (see ``-t`` option of metaphlan3)
              { 'rel_ab', 'rel_ab_w_read_stats', 'reads_map', 'clade_profiles', 'marker_ab_table', 'marker_counts',
              'marker_pres_table', 'clade_specific_strain_tracker' }
        """
        self.analysis_type = self._valid_analysis_type(analysis_type)
        super().__init__(*args, **kwargs)

    def _valid_analysis_type(self, analysis_type):
        choices = [
            "rel_ab", "rel_ab_w_read_stats", "reads_map", "clade_profiles",
            "marker_ab_table", "marker_counts", "marker_pres_table", "clade_specific_strain_tracker"
        ]
        if analysis_type not in choices:
            logger.warning("analysis_type='%s' not valid, set to default ('rel_ab').", analysis_type)
            analysis_type = "rel_ab"
        return analysis_type

    def rows_differences(self, dataframe1, dataframe2) -> pd.DataFrame:
        rows_diff = dataframe1 - dataframe2
        rows_diff[rows_diff.isnull()] = dataframe1
        if self.analysis_type == 'rel_ab':
            rows_diff[rows_diff < 0.0001] = 0
            # if difference between sum of organism of rank r (ex: sum of species of genus X)
            # and value of rank r+1 (ex:genus X) is so small,
            # we assume that it's due to python addition approximation with decimal
        else:
            rows_diff[rows_diff < 0] = 0
        rows_diff = rows_diff.loc[rows_diff.sum(axis=1)[rows_diff.sum(axis=1) != 0].index]
        return rows_diff

    def compare_difference_between_two_levels(self, whole_df, df_at_lower_level, rank) -> pd.DataFrame:
        df_rank = whole_df[whole_df.index.map(lambda x: len(x.split('|'))) == rank]

        # transformation lower_level to rank (level)
        df_rank_computed = df_at_lower_level.copy()
        df_rank_computed.index = df_rank_computed.index.map(lambda x: '|'.join(x.split('|')[:rank]))   # to rank (level)
        df_rank_computed = df_rank_computed.groupby(df_rank_computed.index).sum()             # grouping by rank (level)
        return self.rows_differences(df_rank, df_rank_computed)

    def remove_duplicates(self, df) -> pd.DataFrame:
        """
        Metaphlan3 results are by level therefore we need to remove the duplicated informations
        Example:
        We have:
            ...|g_GenusA    50.0
            ...|g_GenusA|s_Species1 30.0
            ...|g_GenusB    50.0
            ...|g_GenusB|s_Species2 50.0
            Sum = 180.0 =/= 100.0 (while it's relative abundance -> but same problem with other analysis type)
        We want:
            ...|g_GenusA|s_GenusA (genus)   20.0    # unspecified species
            ...|g_GenusA|s_Species1 30.0
            ...|g_GenusB|s_Species2 50.0
            Sum = 100.0
        """
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
                # new_df = new_df.append(rows_diff)              # add missing rows to the dataframe of the lower level
                new_df = pd.concat([new_df, rows_diff])        # add missing rows to the dataframe of the lower level
            # verification that everything is defined up to the lower_level
            samples_with_incomp_lowerlevel = new_df.sum()[new_df.sum() < total]

        new_df = new_df.reset_index()
        return new_df


class Metaphlan2Parser(BaseMetaphlanParser):
    """
    Parse output from `Metaphlan2 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'ID'

    def _load_data(self) -> pd.DataFrame:
        df = super()._load_data()
        df = self.remove_duplicates(df)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df


class Metaphlan3Parser(BaseMetaphlanParser):
    """
    Parse output from `Metaphlan3 <https://github.com/biobakery/MetaPhlAn/>`_ merged table.
    """

    taxa_column = 'clade_name'
    NCBI_tax_column = 'NCBI_tax_id'

    def __init__(self, *args, analysis_type: str = 'rel_ab', keep_NCBI_tax_col: bool = False, **kwargs):
        """
        Args:
            analysis_type: output type of Metaphlan3 (see ``-t`` option of metaphlan3)
              { 'rel_ab', 'rel_ab_w_read_stats', 'reads_map', 'clade_profiles', 'marker_ab_table', 'marker_counts',
              'marker_pres_table', 'clade_specific_strain_tracker' }
            keep_NCBI_tax_col: set to True if you want the NCBI tax column in the returned dataframe.
        """
        self.keep_NCBI_tax_col = keep_NCBI_tax_col
        super().__init__(*args, analysis_type=analysis_type, parsing_options={'skiprows': 1}, **kwargs)

    def _load_data(self) -> pd.DataFrame:
        df = super()._load_data()

        # if number of taxonomical_names is inferior to the default,
        if len(self.taxonomical_names) < len(BaseTaxonomyCountsParser.taxonomical_names):
            # we need to restrict the rows considered to only the rows that recount taxonomical level inside the range
            # wanted.
            # Or error "ValueError: Error : expecting a integer inferior or equal to the number of taxonomical_names."
            # will be raised
            df = df[df["NCBI_tax_id"].map(lambda x: len(x.split("|"))) <= len(self.taxonomical_names)]
        if self.keep_NCBI_tax_col:
            tmp = df[[self.NCBI_tax_column, self.taxa_column]]

        df = df.drop(self.NCBI_tax_column, axis=1)  # NCBI_tax_column needs to be dropped because sum
        df = self.remove_duplicates(df)

        if self.keep_NCBI_tax_col:
            tmp[self.NCBI_tax_column] = tmp[self.NCBI_tax_column].map(lambda x: x.split("|")[-1])
            df = df.merge(tmp)

        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df
