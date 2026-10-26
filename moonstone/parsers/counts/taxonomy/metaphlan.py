import logging
import os
from itertools import takewhile
from typing import List, Optional, Union

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

        Example::

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

    def __init__(
        self, *args, analysis_type: str = 'rel_ab', keep_NCBI_tax_col: bool = False, skiprows: int = 1, **kwargs
    ):
        """
        Args:
            analysis_type: output type of Metaphlan3 (see ``-t`` option of metaphlan3)
              { 'rel_ab', 'rel_ab_w_read_stats', 'reads_map', 'clade_profiles', 'marker_ab_table', 'marker_counts',
              'marker_pres_table', 'clade_specific_strain_tracker' }
            keep_NCBI_tax_col: set to True if you want the NCBI tax column in the returned dataframe.
        """
        self.keep_NCBI_tax_col = keep_NCBI_tax_col
        super().__init__(*args, analysis_type=analysis_type, parsing_options={'skiprows': skiprows}, **kwargs)

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


class MetaphlanReadStatsParser(BaseMetaphlanParser):
    """
    Parse per-sample `MetaPhlAn <https://github.com/biobakery/MetaPhlAn/>`_ profiles made with
    ``-t rel_ab_w_read_stats`` into a table of estimated read counts per clade.

    MetaPhlAn's ``merge_metaphlan_tables.py`` merges only relative abundances, so this parser reads the per-sample
    profiles directly and joins them on clade name. A clade absent from a sample's profile has 0 estimated reads
    in that sample, as in ``merge_metaphlan_tables.py``. The ``UNCLASSIFIED`` row written by MetaPhlAn (reads not
    assigned to a known clade) is kept as its own kingdom.
    """

    taxa_column = 'clade_name'
    reads_columns = ['estimated_number_of_reads_from_the_clade', 'estimated_number_of_bases_from_the_clade']

    def __init__(
        self, file_paths: Union[str, os.PathLike, List[Union[str, os.PathLike]]],
        sample_names: Optional[List[str]] = None, **kwargs
    ):
        """
        Args:
            file_paths: path of one profile, or list of paths (one profile per sample).
            sample_names: names of the samples, in the same order as ``file_paths``. By default, the file name
              without extension and without a ``_profile`` suffix, as in ``merge_metaphlan_tables.py``.
        """
        if isinstance(file_paths, (str, os.PathLike)):
            file_paths = [file_paths]
        self.file_paths = [os.fspath(path) for path in file_paths]
        if sample_names is None:
            sample_names = [
                os.path.splitext(os.path.basename(path))[0].replace('_profile', '') for path in self.file_paths
            ]
        if len(sample_names) != len(self.file_paths):
            raise ValueError(
                f"{len(sample_names)} sample_names given for {len(self.file_paths)} files: expecting one per file."
            )
        duplicates = sorted({name for name in sample_names if list(sample_names).count(name) > 1})
        if duplicates:
            raise ValueError(f"Duplicated sample names: {duplicates}. Give unique sample_names.")
        self.sample_names = list(sample_names)
        super().__init__(self.file_paths, analysis_type='rel_ab_w_read_stats', **kwargs)

    def _read_profile(self, file_path: str) -> pd.Series:
        with open(file_path) as profile:
            headers = list(takewhile(lambda line: line.startswith('#'), profile))
        columns = headers[-1].lstrip('#').rstrip('\n').split('\t') if headers else []
        reads_column = next((column for column in self.reads_columns if column in columns), None)
        if self.taxa_column not in columns or reads_column is None:
            raise ValueError(
                f"{file_path} is not a MetaPhlAn rel_ab_w_read_stats profile: expecting a header line with columns "
                f"'{self.taxa_column}' and '{self.reads_columns[0]}'."
            )
        profile_df = pd.read_csv(
            file_path, sep='\t', skiprows=len(headers), header=None, names=columns,
            usecols=[self.taxa_column, reads_column]
        )
        if profile_df.empty:
            raise ValueError(f"{file_path} has no clade rows.")
        return profile_df.set_index(self.taxa_column)[reads_column]

    def _load_data(self) -> pd.DataFrame:
        df = pd.concat(
            [self._read_profile(path).rename(name) for path, name in zip(self.file_paths, self.sample_names)], axis=1
        ).fillna(0)
        df.index.name = self.taxa_column
        df = df.reset_index()
        # keep only the rows with no more taxonomical levels than taxonomical_names
        df = df[df[self.taxa_column].str.count(r'\|') < len(self.taxonomical_names)]
        df = self.remove_duplicates(df)
        df = self.split_taxa_fill_none(df, sep="|")
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        if (df % 1 == 0).all().all():
            df = df.astype('int64')  # read counts stay integers (pandas turns them to floats while merging)
        return df
