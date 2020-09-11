import pandas as pd
from typing import Optional

from moonstone.plot.graphs.bargraph import (
    BarGraph
)
from moonstone.plot.graphs.histogram import Histogram
from moonstone.utils.plot import (
    add_x_to_plotting_options,
    add_default_titles_to_plotting_options
)


class PlotMetadataStats:
    """
    Several plots available to visualize metadata.
    """

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    def plot_age(self, bins_size=None, plotting_options: dict = None,
                 show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the age distribution of patients (whose the samples are originated from)

        :param bins_size: size of the bins of the Histogram
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {'layout': {'title_text': "Age distribution in the samples", 'title_x': 0.5},
                                'xaxes': {'title_text': "age"},
                                'yaxes': {'title_text': "number of samples"}}
        else:
            plotting_options = add_default_titles_to_plotting_options(plotting_options,
                                                                      "Age distribution in the samples",
                                                                      "age", "number of samples")

        hist_fig = Histogram(self.metadata_df['age'])

        if bins_size is None:
            bins_size = 1

        hist_fig.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file
            )

    def plot_sex(
        self, sex_col: str = "sex", plotting_options: dict = None, show: Optional[bool] = True,
        output_file: Optional[str] = False
    ):
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
            plotting_options = {'layout': {'title_text': "Sex distribution in the samples", 'title_x': 0.5},
                                'xaxes': {'title_text': "sex"},
                                'yaxes': {'title_text': "number of samples"},
                                'traces': {'marker_color': ['pink', 'blue']}}
        else:
            plotting_options = add_default_titles_to_plotting_options(plotting_options,
                                                                      "Sex distribution in the samples",
                                                                      "sex", "number of samples")
            plotting_options = add_x_to_plotting_options(plotting_options, 'traces', 'marker_color', ['pink', 'blue'])

        bar_fig = BarGraph(
            pd.value_counts(self.metadata_df[sex_col]),
        )
        bar_fig.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
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
            plotting_options = {'layout': {'title_text': column_name+" distribution in the samples", 'title_x': 0.5},
                                'xaxes': {'title_text': column_name},
                                'yaxes': {'title_text': "number of samples"}}
        else:
            plotting_options = add_default_titles_to_plotting_options(plotting_options,
                                                                      column_name+" distribution in the samples",
                                                                      column_name, "number of samples")

        bar_fig = BarGraph(
            pd.count_values(self.metadata_df[column_name]),
        )
        # if reset_xnames_dic is not None:
        #     bar_fig.reset_xnames(reset_xnames_dic)
        bar_fig.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            )
