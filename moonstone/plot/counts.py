import pandas as pd
from typing import Optional

from scipy.cluster import hierarchy
from scipy.spatial import distance

from moonstone.filtering import TaxonomyMeanFiltering
from moonstone.plot.graphs.bargraph import BarGraph, MatrixBarGraph
from moonstone.utils.dict_operations import merge_dict
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

    def _compute_top_n_taxa_df(
        self, data_df: pd.DataFrame, taxa_number: int, taxa_level: str
    ):
        df = data_df.groupby(taxa_level).sum()
        top = (
            df[~df.index.str.contains("(", regex=False)]
            .sum(axis=1)
            .sort_values(ascending=False)
            .head(taxa_number)
        )  # Filter out rows not classified to the species level (that contains '(')
        top_and_other_df = df[df.index.get_level_values(taxa_level).isin(top.index)]
        top_and_other_df = df.loc[
            top.index
        ]  # put top species in order from most abundant species across samples to least
        top_and_other_df.loc["Others"] = 100 - top_and_other_df.sum()
        top_and_other_df = top_and_other_df.transpose()

        # Determine samples order using hierarchical clustering
        Z = hierarchy.linkage(
            distance.pdist(
                # top_and_other_df[top_and_other_df.columns[(top_and_other_df>40).any()]].drop('Others', axis=1)),
                # top_and_other_df.drop('Others', axis=1)),
                # top_and_other_df),
                df.transpose()
            ),
            method="single",
            metric="euclidean",
            optimal_ordering=False,
        )
        order = hierarchy.leaves_list(Z)
        return top_and_other_df.iloc[order]

    def plot_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        **kwargs,
    ):
        """
        Plot bar chart of most abundant taxa (total sum of abundance).

        The plot represents percentage of sample with the corresponding taxa
        ordered from most abundant to less abundant.

        Args:
            mean_taxa: mean threshold to be kept for analysis
            taxa_number: number of taxa to plot
            taxa_level: Taxonomy level
        """
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        # Select top `taxa_number`
        top_taxa_mean = (
            data_df.groupby(taxa_level)
            .sum()
            .mean(axis=1)
            .sort_values(ascending=False)[:taxa_number]
        )
        make_float_legend = lambda x: " (mean={:,.2f})".format(x)  # noqa
        top_taxa_mean = top_taxa_mean.apply(make_float_legend)
        number_of_taxa = len(
            top_taxa_mean
        )  # Can be different than threshold if not enough taxa
        # Filter for top species
        abundances = data_df.groupby(taxa_level).sum().loc[top_taxa_mean.index]
        percentage_presence = (abundances != 0).sum(axis=1) / abundances.shape[1] * 100
        percentage_presence.index = percentage_presence.index + top_taxa_mean.astype(
            "str"
        )
        # Make graph
        graph = BarGraph(percentage_presence.iloc[::-1])
        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{number_of_taxa} most abundant {taxa_level} - Total sum of abundances",
                "xaxis_title": "Percentage Sample",
                "yaxis_title": "Species",
            }
        }
        if mean_taxa is not None:
            default_plotting_options["layout"][
                "title"
            ] = "{} (mean among samples > {})".format(
                default_plotting_options["layout"]["title"], mean_taxa
            )
        plotting_options = merge_dict(
            kwargs.pop("plotting_options", {}), default_plotting_options
        )
        graph.plot_one_graph(
            orientation="h", plotting_options=plotting_options, **kwargs
        )

    def plot_sample_composition_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        **kwargs,
    ):
        """
        Plot taxa composition of samples for most abundant taxa.

        Args:
            mean_taxa: mean threshold to be kept for analysis
            taxa_number: number of taxa to plot
            taxa_level: Taxonomy level
        """
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        data_df = self._compute_top_n_taxa_df(data_df, taxa_number, taxa_level)
        df = data_df.T
        # Make graph
        graph = MatrixBarGraph(df)
        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{taxa_level.capitalize()} composition for the top {taxa_number} most abundant species across samples",  # noqa
                "xaxis_title": "Samples",
                "yaxis_title": "Percentage",
                "legend": {"traceorder": "normal"},
            }
        }
        if mean_taxa is not None:
            default_plotting_options["layout"][
                "title"
            ] = "{} (mean among samples > {})".format(
                default_plotting_options["layout"]["title"], mean_taxa
            )
        plotting_options = merge_dict(
            kwargs.pop("plotting_options", {}), default_plotting_options
        )
        graph.plot_one_graph(plotting_options=plotting_options, **kwargs)
