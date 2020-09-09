import pandas as pd
from typing import Optional

from moonstone.plot.plot_template import (
    CategoryBarGraph, DistributionBarGraph, Histogram
)
from moonstone.plot.functions import (
    check_types_in_plotting_options,
    add_x_to_plotting_options
)

"""
plot that could be used on any dataframe no matter the module
Or should we put that in the analysis module?
"""


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
            plotting_options = {}
        else:
            plotting_options = check_types_in_plotting_options(plotting_options)

        plotting_options = add_x_to_plotting_options(
            plotting_options, 'log', True)
        plotting_options = add_x_to_plotting_options(
            plotting_options, 'tickangle', -60)

        df_mean = self.df.mean(axis=1)
        bar_fig = DistributionBarGraph(df_mean, plotting_options, show=show, output_file=output_file)

        bar_fig.compute_heterogeneous_bins()
        bar_fig.in_bins_and_count()    # normalize or not?
        bar_fig.plot_one_graph(
            title="Distribution of %s mean" % self.items_name,
            xlabel="mean of the number of reads",
            ylabel="number of samples",
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
#         if plotting_options is None:
#             plotting_options = {}
#         else:
#             plotting_options = check_types_in_plotting_options(plotting_options)

        # mean by sample or nb of samples with ?

    # heatmap -> in analysis?


class PlotMetadataStats():
    """
    Several plots available to visualize metadata.
    """

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    def plot_age(self, step=None, plotting_options: dict = None,
                 show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the age distribution of patients (whose the samples are originated from)

        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = check_types_in_plotting_options(plotting_options)

        hist_fig = Histogram(self.metadata_df['age'], plotting_options, show=show, output_file=output_file)

        if step is None:
            step = 1

        hist_fig.plot_one_graph(
            "Age distribution in the samples",
            "age",
            "number of samples",
            step
            )

    def plot_sex(self, plotting_options: dict = None, show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the sex distribution of patients (whose the samples are originated from)

        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = check_types_in_plotting_options(plotting_options)

        plotting_options = add_x_to_plotting_options(
            plotting_options, 'colorbar', ['pink', 'blue'])

        bar_fig = CategoryBarGraph(self.metadata_df['sex'], plotting_options, show=show, output_file=output_file)
        # bar_fig.count()    # normalize or not?
        # bar_fig.reset_xnames({'F': 'Female', 'M': 'Male'})
        bar_fig.plot_one_graph(
            title="Sex distribution in the samples",
            xlabel="sex",
            ylabel="number of samples",
            )

    def plot_category_distribution(
        self, column_name, title=None, xlabel=None,
        reset_xnames_dic: dict = None, plotting_options: dict = None,
        show: Optional[bool] = True, output_file: Optional[str] = False
    ):
        """
        :param column_name: name of the column you wish to display into a barplot
        :param title: title of the graph
        :param xlabel: label of the x axis
        :param reset_xnames_dic: to rename the names of the values in the x axis. \n
                                 Example for a plot of the distribution of smoking habits :
                                 reset_xnames_dic={'y': 'smoker', 'n': 'non smoker'}
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = check_types_in_plotting_options(plotting_options)

        bar_fig = CategoryBarGraph(self.metadata_df[column_name], plotting_options, show=show, output_file=output_file)
        bar_fig.count()    # normalize or not?
        if reset_xnames_dic is not None:
            bar_fig.reset_xnames(reset_xnames_dic)
        if title is None:
            title = column_name+" distribution in the samples"
        if xlabel is None:
            xlabel = column_name
        bar_fig.plot_one_graph(
            title,
            xlabel,
            "number of samples",
            )
