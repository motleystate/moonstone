import pandas as pd
from moonstone.utils.pandas.series import SeriesStatsBuilder
from .base import BaseParser
import logging
module_logger = logging.getLogger(__name__)


class MetadataParser(BaseParser):
    """
    Format is the following:

    col_1   col_2   col_3   col_4
    s1  13.3    M   T
    s2  15.3    F   F
    s3  19.1    M   F
    """

    def __init__(self, file_path, sep='\t', no_header=False):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.sep = sep
        self.header = 'infer'
        if no_header is True:
            self.header = None
        super().__init__(file_path)

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path, sep='\t', header=self.header)

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
