import pandas as pd
from typing import Optional, List

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

    def _get_percentage_presence(
        self, df: pd.DataFrame, taxa_level: str, taxa_number: int
    ) -> pd.Series:
        """
        Get percentage presence series for most abundant taxa.

        Args:
            df: Dataframe of abundances
            taxa_level: Taxonomy level
            taxa_number: Number of taxa to plot
        """
        top_taxa_mean = (
            df.groupby(taxa_level)
            .sum()
            .mean(axis=1)
            .sort_values(ascending=False)[:taxa_number]
        )
        make_float_legend = lambda x: " (mean={:,.2f})".format(x)  # noqa
        top_taxa_mean = top_taxa_mean.apply(make_float_legend)
        # Filter for top species
        abundances = df.groupby(taxa_level).sum().loc[top_taxa_mean.index]
        percentage_presence = (abundances != 0).sum(axis=1) / abundances.shape[1] * 100
        percentage_presence.index = percentage_presence.index + top_taxa_mean.astype(
            "str"
        )
        return percentage_presence

    def plot_most_prevalent_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        **kwargs,
    ):
        """
        Plot bar chart of most prevalent taxa (total sum of abundance).

        The plot represents percentage of sample with the corresponding taxa
        ordered from most prevalen to less prevalent.

        Args:
            mean_taxa: Mean threshold to be kept for analysis
            taxa_number: Number of taxa to plot
            taxa_level: Taxonomy level
        """
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        percentage_presence = self._get_percentage_presence(
            data_df, taxa_level, taxa_number
        )
        taxa_number = len(percentage_presence)
        # Make graph
        graph = BarGraph(percentage_presence.iloc[::-1])
        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{taxa_number} most prevalent {taxa_level} - Total sum of abundances",
                "xaxis_title": "Percentage Sample",
                "yaxis_title": taxa_level.capitalize(),
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

    def _compute_top_n_most_abundant_taxa_list(
        self, 
        data_df: pd.DataFrame,
        taxa_number: int
        ):
        return list(
            data_df[~data_df.index.str.contains("(", regex=False)]
            .sum(axis=1)
            .sort_values(ascending=False)
            .head(taxa_number)
            .index
            )  # Filter out rows not classified to the species level (that contains '(')

    def _compute_abundances_taxa_dataframe(
        self, 
        data_df: pd.DataFrame, 
        taxa_level: str,
        taxa_number: int = 20, 
        forced_taxa: list = None,
    ) -> pd.DataFrame:
        """
        Compute abundances for n (taxa_number) taxa with the rest grouped in Others.

        Args:
            data_df: Dataframe of abundances
            taxa_number: Number of taxa to plot
            taxa_level: Taxonomy level
        """
        df = data_df.groupby(taxa_level).sum()

        if forced_taxa:
            top = forced_taxa
        else:
            top = self._compute_top_n_most_abundant_taxa_list(data_df, taxa_number)

        top_and_other_df = df[df.index.get_level_values(taxa_level).isin(top)]
        top_and_other_df = df.loc[
            top
        ]  # put top species in order from most abundant species across samples to least
        top_and_other_df.loc["Others"] = 100 - top_and_other_df.sum()
        return top_and_other_df

    def _cluster_samples(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorder samples with clustering of given order list.

        Args:
            df: dataframe to reorder index
        """
        # Determine samples order using hierarchical clustering
        Z = hierarchy.linkage(
            distance.pdist(df.T),
            method="single",
            metric="euclidean",
            optimal_ordering=False,
        )
        order = hierarchy.leaves_list(Z)
        return df.iloc[:, order]

    def plot_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        cluster_samples: bool = True,
        samples_order: List[str] = None,
        **kwargs,
    ):
        """
        Plot bar chart of most abundant taxa.

        The plot represents percentage of sample with the corresponding taxa
        ordered from most abundant to less abundant.

        Args:
            mean_taxa: Mean threshold to be kept for analysis
            taxa_number: Number of taxa to plot
            taxa_level: Taxonomy level
        """




    def plot_sample_composition_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        cluster_samples: bool = True,
        samples_order: List[str] = None,
        **kwargs,
    ):
        """
        Plot taxa composition of samples for most abundant taxa.

        Args:
            mean_taxa: mean threshold to be kept for analysis
            taxa_number: number of taxa to plot
            taxa_level: Taxonomy level
            cluster_samples: use clustering (skipped by samples_order)
            samples_order: list of samples to force ordering for visualization
        """
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        data_df = self._compute_abundances_for_n_taxa(data_df, taxa_number, taxa_level)
        if samples_order is not None:
            data_df = data_df.loc[:, samples_order]
        elif cluster_samples:
            data_df = self._cluster_samples(data_df)
        # Make graph
        graph = MatrixBarGraph(data_df)
        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{taxa_level.capitalize()} composition for the top {data_df.shape[1]} most abundant species across samples",  # noqa
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
