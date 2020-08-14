import pandas as pd
from typing import Optional

from moonstone.plot.plot_template import (
    Histogram
)

"""
plot that could be used on any dataframe no matter the module
Or should we put that in the analysis module?
"""


# What people might want to visualize?

class PlotStatsData():

    def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
        self.df = dataframe
        self.items_name = items_name

    def _add_x_to_plotting_options_1expectedtype(self, plotting_options: dict, x: str,
                                                 expectedtype: type, defaultvalue):
        if x not in plotting_options.keys():
            plotting_options[x] = defaultvalue
        elif x in plotting_options.keys() and type(plotting_options[x]) != expectedtype:
            raise ValueError('Error : %s value in plotting_options must be a %s, %s given' %
                  (x, expectedtype.__name__,
                   type(plotting_options[x]).__name__))
        return plotting_options

    def _add_x_to_plotting_options_manyexpectedtype(self, plotting_options: dict, x: str,
                                                    expectedtype: list, defaultvalue):
        if x not in plotting_options.keys():
            plotting_options[x] = defaultvalue
        elif x in plotting_options.keys():
            for i in expectedtype:
                if type(plotting_options[x]) == i:
                    return plotting_options
            raise ValueError('Error : %s value in plotting_options must be a %s, %s given' %
                  (x, " or ".join([x.__name__ for x in expectedtype]),
                   type(plotting_options[x]).__name__))
        return plotting_options

    def plot_mean(self, show: Optional[bool] = True, output_file: Optional[str] = False, plotting_options: dict = None):
        """
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options:
        """
        if plotting_options is None:
            plotting_options = {}

        plotting_options = self._add_x_to_plotting_options_1expectedtype(
            plotting_options, 'log', bool, True)
        plotting_options = self._add_x_to_plotting_options_manyexpectedtype(
            plotting_options, 'tickangle', [int, float], -60)

        df_mean = self.df.mean(axis=1)
        hist_fig = Histogram(df_mean, plotting_options, output_file=output_file, show=show)

        hist_fig.compute_asymetric_bins()
        hist_fig.in_bins_and_count()    # normalize or not?
        hist_fig.one_hist(
            "Distribution of %s mean" % self.items_name,
            "mean of the number of reads",
            "number of samples",
            )

    def plot_taxonomy_classification(self, level_of_interest, output_file: Optional[str] = False,
                                     plotting_options: dict = None):
        """
        method to visualize the fractions of taxonomic levels (species/genus/...)
        :param level_of_interest: taxonomic classification level (species, genus etc.) to plot
        :param output_file: name of the output file
        :param plotting_options:
        """
        if plotting_options is None:
            plotting_options = {}

        # mean by sample or nb of samples with ?

    # heatmap -> in analysis?


class PlotStatsMetadata():

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    # plot stats about metadata : proportions M/F etc.?
