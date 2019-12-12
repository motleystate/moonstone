from moonstone.utils.pandas.series import SeriesStatsBuilder


class DataframeStatistics:
    """
    Build statistics on each columns of a given dataframe.
    """

    def __init__(self, dataframe):
        self.dataframe = dataframe

    def get_stats(self):
        """
        :return: list of dict containing statistics about each column
        :rtype: list(dict)
        """
        stats = []
        for col in self.dataframe.columns:
            series_stats_builder = SeriesStatsBuilder(self.dataframe[col])
            stats.append(series_stats_builder.build_stats())
        return stats
