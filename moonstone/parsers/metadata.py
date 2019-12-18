import pandas as pd

from moonstone.analysis.columns_statistics import DataframeStatistics
from .base import BaseParser


class MetadataParser(BaseParser):
    """
    Format is the following:

    col_1   col_2   col_3   col_4
    s1  13.3    M   T
    s2  15.3    F   F
    s3  19.1    M   F
    """

    def __init__(self, file_path, sep='\t', no_header=False):
        self.sep = sep
        self.header = 'infer'
        if no_header is True:
            self.header = None
        super().__init__(file_path)

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path, sep=self.sep, header=self.header)

    def get_stats(self):
        """
        :return: list of dict containing statistics about each column
        :rtype: list(dict)
        """
        return DataframeStatistics(self.dataframe).get_stats()
