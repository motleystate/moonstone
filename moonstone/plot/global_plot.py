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

    def plot_mean(self, output_file: Optional[str] = False, plotting_options: dict = None):
        """
        :param output_file: name of the output file
        :param plotting_options:
        """
        if plotting_options is None:
            plotting_options = {}

        if 'log' in plotting_options.keys() and type(plotting_options['log']) != bool:
            print("Error : log value in plotting_options must be a bool, %s given" %
                  type(plotting_options['log']).__name__)
        elif 'log' not in plotting_options.keys():
            plotting_options['log'] = True

        df_mean = self.df.mean(axis=1)
        hist_fig = Histogram(df_mean, plotting_options, outputfile=output_file)

        binsvalues = hist_fig.compute_asymetric_bins()
        hist_fig.in_bins(binsvalues)
        hist_fig.count()    # normalize or not?
        hist_fig.one_hist(
            "Distribution of %s mean" % self.items_name,
            "mean of the number of reads",
            "number of samples",
            )

    def plot_taxonomy_classification(self, level_of_interest, output_file: Optional[str] = False):
        """
        method to visualize the fractions of taxonomic levels (species/genus/...)
        :param level_of_interest: taxonomic classification level (species, genus etc.) to plot
        :param output_file: name of the output file
        """

    # heatmap -> in analysis?


class PlotStatsMetadata():

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    # plot stats about metadata : proportions M/F etc.?
