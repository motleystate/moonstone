import numpy as np
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
        ordered from most prevalent to less prevalent.

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
            data_df.sum(axis=1)
            .sort_values(ascending=False)
            .head(taxa_number)
            .index
            )

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
        only = False
        if taxa_level[-5:] == "-only":
            taxa_level = taxa_level[:-5]
            only = True

        df = data_df.groupby(taxa_level).sum()

        if forced_taxa:
            top = forced_taxa
        else:
            if only:
                top = self._compute_top_n_most_abundant_taxa_list(
                    df[~df.index.str.contains("(", regex=False)], taxa_number
                    )  # Filter out rows not classified to the species level (that contains '(')
            else:
                top = self._compute_top_n_most_abundant_taxa_list(
                    df, taxa_number
                    )

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

    def _divide_samples_into_subgroups_and_reorder(
        self,
        top_and_other_df,
        sep_series,
        cluster_samples: bool = True
    ):
        ordered_col = []
        x_coor = []
        prec = -0.5
        subgps = sep_series.unique()
        for subgp in subgps:
            if type(subgp) != str and np.isnan(subgp):
                df_gp = top_and_other_df[
                    sep_series[sep_series.isna()].index.intersection(top_and_other_df.columns)]
            else:
                df_gp = top_and_other_df[
                    sep_series[sep_series == subgp].index.intersection(top_and_other_df.columns)]

            if cluster_samples and len(df_gp.columns) > 1:
                tmp = list(self._cluster_samples(df_gp).columns)
                ordered_col += tmp
            else:
                ordered_col += list(df_gp.columns)

            nb = len(df_gp.columns)

            med = nb/2

            x_coor += [(prec, prec+med, prec+nb)]
            # (x of the start of the subgroup annotation square,
            # x of the annotation text,
            # x of the end of the subgroup annotation square)
            prec += nb
        top_and_other_df = top_and_other_df[ordered_col]
        return top_and_other_df, x_coor, subgps

    def _subgroups_annotations(
        self, fig, x_coor, gps
    ):
        i = 0
        color_bg = ["#FFFFFF", "#a7bcdb"]
        while i < len(gps):
            # adding background color
            fig.add_shape(
                type="rect",
                x0=x_coor[i][0], y0=100,
                x1=x_coor[i][2], y1=104,
                line=dict(
                    width=0,
                ),
                fillcolor=color_bg[i % 2],
            )
            # adding text annotation (group name)
            fig.add_annotation(
                x=x_coor[i][1], y=102,
                xref="x", yref="y",
                text=gps[i],
                showarrow=False,
                font=dict(
                    family="Arial",
                    size=14
                    ),
            )
            if i < (len(gps) - 1):
                # adding line separating groups
                fig.add_shape(
                    type="line",
                    x0=x_coor[i][2], y0=100,
                    x1=x_coor[i][2], y1=0,
                    line=dict(
                        width=1,
                        dash="solid",
                        color="white"
                    )
                )
            i += 1
        return fig

    def plot_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
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
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df

        abundance_per_samples_df = self._compute_abundances_taxa_dataframe(data_df, taxa_level, taxa_number=taxa_number)
        average_abundance_df = abundance_per_samples_df.transpose().drop('Others', axis=1).mean()

        taxa_number = len(average_abundance_df)
        # Make graph
        graph = BarGraph(average_abundance_df)
        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{taxa_number} most abundant {taxa_level}",
                "xaxis_title": "Average relative abundance",
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

    def plot_sample_composition_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        cluster_samples: bool = True,
        samples_order: List[str] = None,
        color_df: pd.DataFrame = None,
        sep_series: pd.Series = None,
        sep_how: str = None,
        **kwargs,
    ):
        """
        Plot taxa composition of samples for most abundant taxa.

        Args:
            mean_taxa: mean threshold to be kept for analysis
            taxa_number: number of taxa to plot
            taxa_level: Taxonomy level. Add "-only" after the taxonomy level if you do not want OTU only defined
            at a higher level to appear in the top. They will still appear in "Others"
            cluster_samples: use clustering (skipped by samples_order)
            samples_order: list of samples to force ordering for visualization
            color_df: metadata to put as legend on the bottom of the graph
            sep_series: metadata used to order samples into subgroups (skipped by samples_order)
            sep_how: { None, 'color', 'labels' } Graphical way of showing the separation of the different subgroups
            (skipped if sep_series is empty/None)
        """
        data_df = self.df
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        data_df = self._compute_abundances_taxa_dataframe(data_df, taxa_level, taxa_number=taxa_number)

        if data_df.shape[1] <= 1:        # only 1 sample, no need for ordering
            sep_series = None
        elif samples_order is not None:
            data_df = data_df.loc[:, samples_order]
        elif sep_series is not None:     # organize samples inside subgroups and concatenate subgroups one after another
            data_df, x_coor, subgps = self._divide_samples_into_subgroups_and_reorder(
                data_df, sep_series, cluster_samples=cluster_samples
                )
            if sep_how == 'color':
                if color_df is None:
                    color_df = pd.DataFrame(sep_series)
                elif sep_series.name not in color_df.columns:
                    color_df = color_df.merge(pd.DataFrame(sep_series), right_index=True, left_index=True)
        elif cluster_samples:
            data_df = self._cluster_samples(data_df)

        # Make graph
        graph = MatrixBarGraph(data_df)
        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{taxa_level.capitalize()} composition for the top {data_df.drop('Others').shape[0]} most abundant species across samples",  # noqa
                "xaxis_title": "Samples",
                "yaxis_title": "Percentage",
                "legend": {"traceorder": "normal"},
                "legend_title_text": "species",
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

        if sep_series is not None and sep_how == "labels":
            show = kwargs.pop("show", True)
            output_file = kwargs.pop("output_file", False)
            if color_df is None:
                fig = graph.plot_one_graph(plotting_options=plotting_options, **kwargs, show=False)
            else:
                fig = graph.plot_complex_graph(color_df, plotting_options=plotting_options, **kwargs, show=False)
            fig = self._subgroups_annotations(fig, x_coor, subgps)
            graph._handle_output_plotly(fig, show, output_file)
        else:
            if color_df is None:
                fig = graph.plot_one_graph(plotting_options=plotting_options, **kwargs)
            else:
                fig = graph.plot_complex_graph(color_df, plotting_options=plotting_options, **kwargs)
