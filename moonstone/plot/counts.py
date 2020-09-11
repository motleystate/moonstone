import pandas as pd
from typing import Optional

from moonstone.plot.graphs.bargraph import BarGraph
from moonstone.utils.plot import (
    add_x_to_plotting_options,
    add_default_titles_to_plotting_options
)
from moonstone.utils.pandas.series import SeriesBinning


# What people might want to visualize?


class PlotCountsStats():
    """
    Several plots available to visualize simple count data.
    """

    def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
        self.df = dataframe
        self.items_name = items_name

    def plot_mean_distribution(
        self, plotting_options: dict = None, show: Optional[bool] = True, output_file: Optional[str] = False
    ):
        """
        method to visualize the mean distribution of the number of reads by items

        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {'layout': {'title_text': "Distribution of %s mean" % self.items_name, 'title_x': 0.5,
                                           'yaxis_type': 'log'},
                                'xaxes': {'title_text': "sex", 'tickangle': -60},
                                'yaxes': {'title_text': "number of samples"}}
        else:
            plotting_options = add_default_titles_to_plotting_options(plotting_options,
                                                                      "Distribution of %s mean" % self.items_name,
                                                                      "mean of the number of reads",
                                                                      "number of samples")
            plotting_options = add_x_to_plotting_options(plotting_options, 'layout', 'yaxis_type', 'log')
            plotting_options = add_x_to_plotting_options(plotting_options, 'xaxes', 'tickangle', -60)

        mean_series = self.df.mean(axis=1)
        binned_mean = SeriesBinning(mean_series).binned_data
        bar_fig = BarGraph(binned_mean, plotting_options, show=show, output_file=output_file)
        bar_fig.plot_one_graph(
            plotting_options,
            show=show,
            output_file=output_file
        )


# class PlotTaxonomyStats():
#     """
#     Several plots available to visualize taxonomy count data.
#     """

#     def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
#         self.df = dataframe
#         self.items_name = items_name

#     def _plot_taxonomy_classification(self, level_of_interest, plotting_options: dict = None,
#                                       show: Optional[bool] = True, output_file: Optional[str] = False):
#         """
#         ~~~~ IN CONSTRUCTION (remove _ in title when finished) ~~~~
#         method to visualize the fractions of taxonomic levels (species/genus/...)

#         :param level_of_interest: taxonomic classification level (species, genus etc.) to plot
#         :param output_file: name of the output file
#         :param plotting_options: options of plotting that will override the default setup \n
#                                  [!] Make sure the value given to an argument is of the right type \n
#                                  options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
#                                  'tickangle': `[int, float]`
#         """

        # mean by sample or nb of samples with ?

    # heatmap -> in analysis?
