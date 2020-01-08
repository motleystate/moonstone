from moonstone.utils.convert import pandas_to_python_type


class SeriesStatsBuilder(object):

    def __init__(self, series):
        self.series = series

    def _build_base_stats(self):
        repartition = self.series.value_counts()
        return {
            'col_name': self.series.name,
            'col_type': str(self.series.dtype),
            'python_col_type': pandas_to_python_type(self.series.dtype),
            'n_values': self.series.size,
            'n_uniq_values': len(self.series.value_counts()),
            'values_repartition': {i: repartition[i] for i in repartition.index},
        }

    def _build_object_stats(self):
        return self._build_base_stats()

    def _build_number_stats(self):
        stats_dict = self._build_base_stats()
        stats_dict['mean'] = round(self.series.mean(), 2)
        return stats_dict

    def _build_int64_stats(self):
        return self._build_number_stats()

    def _build_float64_stats(self):
        return self._build_number_stats()

    def build_stats(self):
        return getattr(self, f"_build_{str(self.series.dtype)}_stats",
                       self._build_base_stats)()
