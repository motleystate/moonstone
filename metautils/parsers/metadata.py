import pandas as pd

from metautils.utils.convert import pandas_to_python_type

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

    def _build_base_stats(self, series):
        return {
            'col_name': series.name,
            'col_type': str(series.dtype),
            'python_col_type': pandas_to_python_type(series.dtype),
            'n_values': series.size,
            'n_uniq_values': len(series.value_counts())
        }

    def _build_object_stats(self, series):
        stats_dict = self._build_base_stats(series)
        repartition = series.value_counts()
        stats_dict['values_repartition'] = {i: repartition[i] for i in repartition.index}
        return stats_dict

    def _build_number_stats(self, series):
        stats_dict = self._build_base_stats(series)
        stats_dict['mean'] = round(series.mean(), 2)
        return stats_dict

    def _build_int64_stats(self, series):
        return self._build_number_stats(series)

    def _build_float64_stats(self, series):
        return self._build_number_stats(series)

    def _build_stats(self, col_name):
        """
        Build statistics of a pandas Series from dataframe based on column name
        :return: Dict of stats
        :rtype: dict
        """
        series = self.dataframe[col_name]
        return getattr(self, f"_build_{str(series.dtype)}_stats",
                       self._build_base_stats)(series)

    def get_stats(self):
        """
        :return: list of dict containing statistics about each column
        :rtype: list(dict)
        """
        stats = []
        for i in self.dataframe.columns:
            stats.append(self._build_stats(i))
        return stats
