import pandas as pd
from typing import Optional

from moonstone.filtering import TaxonomyMeanFiltering
from moonstone.plot.graphs.bargraph import BarGraph
from moonstone.utils.plot import (
    add_x_to_plotting_options,
    add_default_titles_to_plotting_options,
)
from moonstone.utils.pandas.series import SeriesBinning


# What people might want to visualize?


class PlotCountsStats:
    """
    Several plots available to visualize simple count data.
    """

    def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
        self.df = dataframe
        self.items_name = items_name

    def plot_mean_distribution(
        self,
        plotting_options: dict = None,
        show: Optional[bool] = True,
        output_file: Optional[str] = False,
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
            plotting_options = {
                "layout": {
                    "title_text": "Distribution of %s mean" % self.items_name,
                    "title_x": 0.5,
                    "yaxis_type": "log",
                },
                "xaxes": {"title_text": "sex", "tickangle": -60},
                "yaxes": {"title_text": "number of samples"},
            }
        else:
            plotting_options = add_default_titles_to_plotting_options(
                plotting_options,
                "Distribution of %s mean" % self.items_name,
                "mean of the number of reads",
                "number of samples",
            )
            plotting_options = add_x_to_plotting_options(
                plotting_options, "layout", "yaxis_type", "log"
            )
            plotting_options = add_x_to_plotting_options(
                plotting_options, "xaxes", "tickangle", -60
            )

        mean_series = self.df.mean(axis=1)
        binned_mean = SeriesBinning(mean_series).binned_data
        bar_fig = BarGraph(
            binned_mean, plotting_options, show=show, output_file=output_file
        )
        bar_fig.plot_one_graph(plotting_options, show=show, output_file=output_file)


class PlotTaxonomyCounts:
    """
    Plots available for taxonomy counts (multiindexed dataframe).
    """

    def __init__(self, taxonomy_dataframe: pd.DataFrame):
        self.df = taxonomy_dataframe

    def plot_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        **kwargs,
    ):
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        # Select top `taxa_number`
        top_taxa = (
            data_df.groupby(taxa_level)
            .sum()
            .sum(axis=1)
            .sort_values(ascending=False)[:taxa_number]
            .index
        )
        number_of_taxa = len(
            top_taxa
        )  # Can be different than threshold if not enough taxa
        # Filter for top species
        abundances = data_df.groupby(taxa_level).sum().loc[top_taxa]
        percentage_presence = (abundances != 0).sum(axis=1) / abundances.shape[1] * 100
        # Make graph
        graph = BarGraph(percentage_presence.iloc[::-1])
        # Plotting options
        plotting_options = {
            "layout": {
                "title": f"{number_of_taxa} most abundant {taxa_level} - Total sum of abundances",
                "xaxis_title": "Percentage Sample",
                "yaxis_title": "Species",
            }
        }
        if mean_taxa is not None:
            plotting_options["layout"]["title"] = "{} (mean among samples > {})".format(
                plotting_options["layout"]["title"], mean_taxa
            )
        graph.plot_one_graph(
            orientation="h", plotting_options=plotting_options, **kwargs
        )
