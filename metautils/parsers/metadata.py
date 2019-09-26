import pandas as pd

from metautils.utils.pandas.series import SeriesStatsBuilder

from .base import BaseParser


class MetadataParser(BaseParser):
    """
    Format is the following:

    col_1   col_2   col_3   col_4
    s1  13.3    M   T
    s2  15.3    F   F
    s3  19.1    M   F
    """

    def __init__(self, file_path, sep='\t'):
        self.sep = sep
        super().__init__(file_path)

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path, sep='\t')

    def get_stats(self):
        """
        :return: list of dict containing statistics about each column
        :rtype: list(dict)
        """
        stats = []
        for col in self.dataframe.columns:
            stats_builder = SeriesStatsBuilder(self.dataframe[col])
            stats.append(stats_builder.build_stats())
        return stats
