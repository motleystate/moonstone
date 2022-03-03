import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from typing import Optional, List, Tuple, Union

from scipy.cluster import hierarchy
from scipy.spatial import distance

from moonstone.plot.graphs.bargraph import BarGraph, MatrixBarGraph
from moonstone.plot.graphs.box import GroupBoxGraph
from moonstone.plot.graphs.violin import GroupViolinGraph
from moonstone.utils.dict_operations import merge_dict
from moonstone.utils.plot import (
    add_x_to_plotting_options,
    add_default_titles_to_plotting_options,
    add_groups_annotations,
)
from moonstone.utils.pandas.series import SeriesBinning

logger = logging.getLogger(__name__)

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
            plotting_options = add_x_to_plotting_options(plotting_options, "layout", "yaxis_type", "log")
            plotting_options = add_x_to_plotting_options(plotting_options, "xaxes", "tickangle", -60)

        mean_series = self.df.mean(axis=1)
        binned_mean = SeriesBinning(mean_series).binned_data
        bar_fig = BarGraph(binned_mean, plotting_options, show=show, output_file=output_file)
        bar_fig.plot_one_graph(plotting_options, show=show, output_file=output_file)


class PlotTaxonomyCounts:
    """
    Plots available for taxonomy counts (multiindexed dataframe).
    """

    def __init__(self, taxonomy_dataframe: pd.DataFrame):
        if isinstance(taxonomy_dataframe, pd.Series):
            self.df = pd.DataFrame(taxonomy_dataframe)
        else:
            self.df = taxonomy_dataframe

    def compute_prevalence_series(self) -> pd.Series:
        return (self.df != 0).sum(axis=1) / self.df.shape[1] * 100

    @property
    def prevalence_series(self):
        # call compute_prevalence_series and store into self._prevalence_series
        if getattr(self, "_prevalence_series", None) is None:
            self._prevalence_series = self.compute_prevalence_series()
        return self._prevalence_series

    def compute_relative_abundance_dataframe(self) -> pd.DataFrame:
        return self.df * 100 / self.df.sum()

    @property
    def relative_abundance_dataframe(self):
        # call compute_relative_abundance_dataframe and store into self._relative_abundance_dataframe
        if getattr(self, "_relative_abundance_dataframe", None) is None:
            self._relative_abundance_dataframe = self.compute_relative_abundance_dataframe()
        return self._relative_abundance_dataframe

    def _valid_mode_param(self, mode: str) -> str:
        if mode[:3] == "box":
            return "boxplot"
        if mode[:6] == "violin":
            return "violin"
        if mode[:3] == "bar":
            return "bargraph"
        logger.warning("mode='%s' not valid, set to default (bargraph).", mode)
        return "bargraph"

    def _add_mean_info_to_index(self, data_df: Union[pd.DataFrame, pd.Series], mean_counts_ser_taxa: pd.Series):
        """
        add the mean of the taxa among samples at the end of its name in the index.
        """
        mean_top_prev = mean_counts_ser_taxa.loc[data_df.index]
        make_float_legend = lambda x: " (mean={:,.2f})".format(x)  # noqa
        mean_top_prev = mean_top_prev.apply(make_float_legend)
        data_df.index = data_df.index + mean_top_prev.astype("str")
        return data_df

    def _italicize_taxa_name(self, text: str) -> str:
        """
        put <i> and </i> around the taxa name, so that it will be shown in italic in the graph.
        It leaves the mean info or the higher classification between parenthesis, unitalicized.

        Args:
            text: string that contains taxa name to italicize
        """
        s = text.split(" (")
        taxa_name = s[0].replace("_", " ")
        taxa_name = "<i>" + taxa_name + "</i>"
        if len(s) > 1:
            end = " (" + " (".join(s[1:])
            return taxa_name + end
        else:
            return taxa_name

    def _generate_list_species_to_plot(
        self,
        determining_ser_taxa: pd.Series,
        other_variable_ser_taxa: pd.Series,
        taxa_number: int = 20,
        determining_threshold: float = None,
        higher_classification: bool = True,
        threshold_on_other_variable: float = None,
        ascending: bool = False,
    ) -> pd.core.indexes.base.Index:
        """
        generate the list of species to plot, the most abundant/prevalent

        Args:
            determining_ser_taxa: Series used to compute the top most abundant/prevalent taxa.
            other_variable_ser_taxa: Series used to filter some taxa out based on the other statistical variable.
              (so on prevalence, if plotting most abundant taxa, and on mean counts, if plotting most prevalent taxa)
            taxa_number: Number of taxa to plot (skipped by determining_threshold).
            determining_threshold: (optional) Set a threshold, if rather than show a certain number of taxa, you want
              to show all taxa with an equal or higher prevalence/relative abundance.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top. They will still be included in relative abundances.
            threshold_on_other_variable: (optional) The threshold of the other variable used to filter some taxa out.
            ascending: Set to True, if you want the taxa to be ordered from least abundant/prevalent taxa of the top,
              to most abundant/prevalent taxa of the top.
        """
        # determining_ser_taxa: mean_relab_ser_taxa for abundance; prev_ser_taxa for prevalence
        # other_variable_ser_taxa: prev_ser_taxa for abundance; mean_counts_ser_taxa for prevalence
        # determining_threshold: average_relative_abundance_threshold for abundance; prevalence_threshold for prevalence
        # threshold_on_other_variable: prevalence_threshold for abundance; mean_threshold for prevalence

        if not higher_classification:
            other_variable_ser_taxa = other_variable_ser_taxa[
                ~other_variable_ser_taxa.index.str.contains("(", regex=False)
            ]
        if threshold_on_other_variable:
            sp_to_keep = other_variable_ser_taxa[other_variable_ser_taxa >= threshold_on_other_variable].index
        else:
            sp_to_keep = other_variable_ser_taxa.index

        determining_ser_taxa = determining_ser_taxa.loc[sp_to_keep]
        if determining_threshold:
            top_species = determining_ser_taxa[determining_ser_taxa >= determining_threshold]
            top_species = top_species.sort_values(ascending=ascending).index
        else:
            top_species = determining_ser_taxa.sort_values(ascending=False)[:taxa_number].index
            if ascending:
                top_species = top_species[::-1]

        if len(top_species) == 0:
            logger.warning(
                "No species abide by the threshold(s) given. You may want to try to lower your threshold(s)."
            )
        return top_species

    def _cluster_samples(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorder samples with clustering of given order list.

        Args:
            df: dataframe whose index needs to be reordered
        """
        # Determine samples order using hierarchical clustering
        Z = hierarchy.linkage(
            distance.pdist(df.drop(["Others"]).T),
            method="single",
            metric="euclidean",
            optimal_ordering=False,
        )
        order = hierarchy.leaves_list(Z)
        return df.iloc[:, order]

    def _divide_samples_into_subgroups_and_reorder(
        self,
        top_and_other_df: pd.DataFrame,
        sep_series: pd.Series,
        cluster_samples: bool = True,
    ) -> Tuple[pd.DataFrame, List[Tuple[float, float, float]], np.ndarray]:
        """
        divide samples into subgroups and reorder them using hierarchical clustering if cluster_samples is True.

        Args:
            top_and_other_df: DataFrame of the relative abundances of the top taxa and of Others (every taxa not in the
              top is considered as Others).
            sep_series: Series of the metadata used to divide samples.
            cluster_samples: Set to False, if you don't want the samples to be clusterize using hierarchical clustering.
        """
        ordered_col = []
        x_coor = []
        prec = -0.5
        subgps = sep_series.unique()
        for subgp in subgps:
            if type(subgp) != str and np.isnan(subgp):
                df_gp = top_and_other_df[sep_series[sep_series.isna()].index.intersection(top_and_other_df.columns)]
            else:
                df_gp = top_and_other_df[sep_series[sep_series == subgp].index.intersection(top_and_other_df.columns)]

            if cluster_samples and len(df_gp.columns) > 1:
                tmp = list(self._cluster_samples(df_gp).columns)
                ordered_col += tmp
            else:
                ordered_col += list(df_gp.columns)

            nb = len(df_gp.columns)

            med = nb / 2

            x_coor += [(prec, prec + med, prec + nb)]
            # (x of the start of the subgroup annotation square,
            # x of the annotation text,
            # x of the end of the subgroup annotation square)
            prec += nb
        top_and_other_df = top_and_other_df[ordered_col]
        return top_and_other_df, x_coor, subgps

    def _compute_relative_abundances_taxa_dataframe(
        self,
        taxa_level: str = "species",
        taxa_number: int = 20,
        average_relative_abundance_threshold: float = None,
        higher_classification: bool = True,
        prevalence_threshold: float = None,
        ascending: bool = False,
    ) -> pd.DataFrame:
        """
        Compute for each samples, the relative abundances for the n (taxa_number) most abundant taxa across all the
        samples.

        Args:
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot (skipped by average_relative_abundance_threshold).
            average_relative_abundance_threshold: (optional) Set a threshold, if you want to show all taxa with an
              equal or greater average relative abundance.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top. They will still be included in the relative abundances.
            prevalence_threshold: Prevalence threshold for a taxa to be kept in analysis.
            ascending: If set to True, from top to bottom, from least prevalent taxa of the top to most prevalent taxa.
        """
        relab_df_taxa = self.relative_abundance_dataframe.groupby(taxa_level).sum()
        # if taxa isn't the lowest taxonomical level, it sums up all counts of the same taxa
        # if taxa is the lowest taxonomical level, it drops the higher taxonomical levels in index
        # MultiIndex -> (single) Index

        mean_relab_ser_taxa = relab_df_taxa.mean(axis=1)

        prev_ser_taxa = self.prevalence_series.groupby(taxa_level).mean()

        top_ab = self._generate_list_species_to_plot(
            mean_relab_ser_taxa,
            prev_ser_taxa,
            taxa_number,
            average_relative_abundance_threshold,
            higher_classification,
            prevalence_threshold,
            ascending,
        )

        taxa_number = len(top_ab)  # for prevalence_thresholds case
        # and also in the case that there is less species that the number asked for

        return relab_df_taxa.loc[top_ab], taxa_number, mean_relab_ser_taxa.loc[top_ab]

    def _plot_most_abundant_taxa_bargraph(
        self,
        taxa_level: str = "species",
        taxa_number: int = 20,
        average_relative_abundance_threshold: float = None,
        higher_classification: bool = True,
        prevalence_threshold: float = None,
        ascending: bool = False,
        plotting_options: dict = {},
        **kwargs,
    ) -> go.Figure:
        """
        Generate Bar Graph of the most abundant taxa

        Args:
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot (skipped by average_relative_abundance_threshold).
            average_relative_abundance_threshold: (optional) Set a threshold, if you want to show all species with an
              equal or greater average relative abundance.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top. They will still be included in relative abundances.
            prevalence_threshold: Prevalence threshold for a taxa to be kept in analysis.
            ascending: If set to True, from top to bottom, from least abundant taxa of the top to most abundant taxa.
            plotting_options: options for the layout of the graph.
        """
        (taxa_number, average_abundance_ser,) = self._compute_relative_abundances_taxa_dataframe(
            taxa_level=taxa_level,
            taxa_number=taxa_number,
            average_relative_abundance_threshold=average_relative_abundance_threshold,
            higher_classification=higher_classification,
            prevalence_threshold=prevalence_threshold,
            ascending=bool(1 - ascending),
        )[1:]
        average_abundance_ser.index = average_abundance_ser.index.map(self._italicize_taxa_name)

        # Make graph
        graph = BarGraph(average_abundance_ser)
        # Plotting options
        title = f"{taxa_number} most abundant {taxa_level}"
        if prevalence_threshold:
            title += f" (present in at least {prevalence_threshold}% of samples)"
        default_plotting_options = {
            "layout": {
                "title": title,
                "xaxis_title": "Average relative abundance",
                "yaxis_title": taxa_level.capitalize(),
            }
        }

        plotting_options = merge_dict(plotting_options, default_plotting_options)

        fig = graph.plot_one_graph(orientation="h", plotting_options=plotting_options, **kwargs)

        return fig

    def _plot_most_what_taxa_boxplot_or_violin(
        self,
        what: str,
        mode: str,
        taxa_level: str = "species",
        taxa_number: int = 20,
        determining_threshold: float = None,
        higher_classification: bool = True,
        threshold_on_other_variable: bool = False,
        ascending: bool = False,
        plotting_options: dict = {},
        mean_info: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Generate Box or Violin plot showing every samples' relative abundance as a point, for the top most
        abundant/prevalent taxa.

        Args:
            what: { "abundant", "prevalent" } Variable used to determine and order the taxa shown in the graph.
            mode: { "boxplot", "violin" } Mode of the graph .
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot (skipped by determining_threshold).
            determining_threshold: (optional) Set a threshold, if rather than show a certain number of taxa, you want
              to show all taxa with an equal or higher prevalence/relative abundance.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top. They will still be included in relative abundances.
            threshold_on_other_variable: (optional) Set a threshold to filter some taxa out of top, based on the other
              variable.
            ascending: If set to True, from top to bottom, from least abundant/prevalent taxa of the top to most
              abundant/prevalent taxa.
            mean_info: To show the taxa's mean counts in the graph, next to the taxa's name.
            plotting_options: options for the layout of the graph.
        """
        # The kind of plot for both most abundant and most prevalent, except for the taxa represented

        title = ""
        ascending = bool(1 - ascending)
        # to have from top to bottom, the most abundant/prevalent species to the least abundant/prevalent species of the
        # top, you need to sort_values descendingly

        relab_df_taxa = self.relative_abundance_dataframe.groupby(taxa_level).sum()
        # if taxa isn't the lowest taxonomical level, it sums up all counts of the same taxa
        # if taxa is the lowest taxonomical level, it drops the higher taxonomical levels in index
        # MultiIndex -> (single) Index

        prev_ser_taxa = self.prevalence_series.groupby(taxa_level).mean()

        if what == "abundant":
            mean_relab_ser_taxa = relab_df_taxa.mean(axis=1)

            groups = self._generate_list_species_to_plot(
                mean_relab_ser_taxa,
                prev_ser_taxa,
                taxa_number,
                determining_threshold,  # determining_threshold = average_relative_abundance_threshold
                higher_classification,
                threshold_on_other_variable,  # threshold_on_other_variable = prevalence_threshold
                ascending,
            )

            relab_df_taxa = relab_df_taxa.loc[groups]

            if threshold_on_other_variable:  # threshold_on_other_variable = prevalence_threshold
                title = f" (present in at least {threshold_on_other_variable}% of samples)"

        if what == "prevalent":
            mean_counts_ser_taxa = self.df.groupby(taxa_level).sum().mean(axis=1)

            groups = self._generate_list_species_to_plot(
                prev_ser_taxa,
                mean_counts_ser_taxa,
                taxa_number,
                determining_threshold,  # determining_threshold = prevalence_threshold
                higher_classification,
                threshold_on_other_variable,  # threshold_on_other_variable = mean_threshold
                ascending=ascending,
            )

            relab_df_taxa = relab_df_taxa.loc[groups]

            if mean_info:
                # adding mean information
                relab_df_taxa = self._add_mean_info_to_index(relab_df_taxa, mean_counts_ser_taxa)
                groups = list(relab_df_taxa.index)

            if threshold_on_other_variable:  # threshold_on_other_variable = mean_threshold
                title = f" (with mean among samples > {threshold_on_other_variable})"

        nb = relab_df_taxa.shape[0]
        relab_df_taxa2 = relab_df_taxa[relab_df_taxa.columns[0]].reset_index()
        relab_df_taxa2.index = nb * [relab_df_taxa.columns[0]]
        relab_df_taxa2.columns = ["species", "relative abundance"]
        for i in relab_df_taxa.columns[1:]:
            tmp = relab_df_taxa[i].reset_index()
            tmp.index = nb * [i]
            tmp.columns = ["species", "relative abundance"]
            relab_df_taxa2 = relab_df_taxa2.append(tmp)
        relab_df_taxa2.species = relab_df_taxa2.species.apply(self._italicize_taxa_name)
        groups = [self._italicize_taxa_name(name) for name in groups]

        # Make graph
        if mode == "violin":
            graph = GroupViolinGraph(relab_df_taxa2)
        else:
            graph = GroupBoxGraph(relab_df_taxa2)
        # Plotting options
        final_colors = {name: "#778899" for name in groups}

        default_plotting_options = {
            "layout": {
                "title": f"Relative abundance of the {len(groups)} most {what} microbial genomes among individuals \
of the cohort"
                + title,
                "xaxis_type": "log",
                "showlegend": False,
            },
            "xaxes": {"title_text": "Relative abundance (in percentage)"},
        }

        plotting_options = merge_dict(plotting_options, default_plotting_options)

        fig = graph.plot_one_graph(
            data_col="relative abundance",
            group_col="species",
            groups=groups,
            colors=final_colors,
            orientation="h",
            plotting_options=plotting_options,
            **kwargs,
        )

        return fig

    def _plot_most_prevalent_taxa_bargraph(
        self,
        taxa_level: str = "species",
        taxa_number: int = 20,
        prevalence_threshold: float = None,
        higher_classification: bool = True,
        mean_threshold: float = None,
        mean_info: bool = True,
        ascending: bool = False,
        plotting_options: dict = {},
        **kwargs,
    ) -> go.Figure:
        """
        Generate Bar Graph of the most prevalent taxa

        Args:
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot (skipped by prevalence_threshold).
            prevalence_threshold: (optional) Set a threshold, if you want to show all species with an equal or greater
              prevalence.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top.
            mean_threshold: Mean threshold for a taxa to be kept in analysis.
            mean_info: To show the taxa's mean counts in the graph, next to the taxa's name.
            ascending: If set to True, from top to bottom, from least prevalent taxa of the top to most prevalent taxa.
            plotting_options: options for the layout of the graph.
        """
        ascending = bool(1 - ascending)
        title = ""

        prev_ser_taxa = self.prevalence_series.groupby(taxa_level).mean()
        # if taxa isn't the lowest taxonomical level, it sums up all counts of the same taxa
        # if taxa is the lowest taxonomical level, it drops the higher taxonomical levels in index
        # MultiIndex -> (single) Index

        mean_counts_ser_taxa = self.df.groupby(taxa_level).sum().mean(axis=1)

        top_prev = self._generate_list_species_to_plot(
            prev_ser_taxa,
            mean_counts_ser_taxa,
            taxa_number,
            prevalence_threshold,
            higher_classification,
            mean_threshold,
            ascending,
        )

        prev_ser_taxa = prev_ser_taxa.loc[top_prev]

        if mean_info:
            # adding mean information
            prev_ser_taxa = self._add_mean_info_to_index(prev_ser_taxa, mean_counts_ser_taxa)
        prev_ser_taxa.index = prev_ser_taxa.index.map(self._italicize_taxa_name)

        if mean_threshold:
            title = f" (with mean among samples > {mean_threshold})"

        taxa_number = len(top_prev)  # in the case that there is less species that the number asked for

        # Make graph
        graph = BarGraph(prev_ser_taxa)

        # Plotting options
        default_plotting_options = {
            "layout": {
                "title": f"{taxa_number} most prevalent {taxa_level}" + title,
                "xaxis_title": "Percentage Sample",
                "yaxis_title": taxa_level.capitalize(),
            }
        }

        plotting_options = merge_dict(plotting_options, default_plotting_options)

        fig = graph.plot_one_graph(orientation="h", plotting_options=plotting_options, **kwargs)

        return fig

    def plot_most_prevalent_taxa(
        self,
        mode: str = "bargraph",
        taxa_level: str = "species",
        taxa_number: int = 20,
        prevalence_threshold: float = None,
        higher_classification: bool = True,
        mean_threshold: float = None,
        mean_info: bool = False,
        ascending: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Generate a plot of most prevalent taxa.

        Args:
            mode: { 'bargraph' (default), 'boxplot', 'violin' } Bargraph will show you the prevalence of the most
              prevalent taxa among all the samples. Boxplot and violin plot will show every samples' relative abundance
              as a point, for the top most prevalent taxa.
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot (skipped by prevalence_threshold).
            prevalence_threshold: (optional) Set a threshold, if you want to show all species with an equal or greater
              prevalence.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top.
            mean_threshold: Mean threshold for a taxa to be kept in analysis.
            mean_info: To show the taxa's mean counts in the graph, next to the taxa's name.
            ascending: If set to True, from top to bottom, from least prevalent taxa of the top to most prevalent taxa.
        """
        plotting_options = kwargs.pop("plotting_options", {})
        mode = self._valid_mode_param(mode)
        if mode == "bargraph":
            fig = self._plot_most_prevalent_taxa_bargraph(
                taxa_level=taxa_level,
                taxa_number=taxa_number,
                prevalence_threshold=prevalence_threshold,
                higher_classification=higher_classification,
                mean_threshold=mean_threshold,
                mean_info=mean_info,
                ascending=ascending,
                plotting_options=plotting_options,
                **kwargs,
            )
        else:
            fig = self._plot_most_what_taxa_boxplot_or_violin(
                "prevalent",
                mode,
                taxa_level=taxa_level,
                taxa_number=taxa_number,
                determining_threshold=prevalence_threshold,
                higher_classification=higher_classification,
                threshold_on_other_variable=mean_threshold,
                ascending=ascending,
                plotting_options=plotting_options,
                mean_info=mean_info,
                **kwargs,
            )

        return fig

    def plot_most_abundant_taxa(
        self,
        mode: str = "bargraph",
        taxa_level: str = "species",
        taxa_number: int = 20,
        average_relative_abundance_threshold: float = None,
        higher_classification: bool = True,
        prevalence_threshold: float = None,
        ascending: bool = False,
        **kwargs,
    ) -> go.Figure:
        """
        Generate a plot of most abundant taxa.

        Args:
            mode: { 'bargraph' (default), 'boxplot', 'violin' } Bargraph will show you the mean relative abundance of
              the most abundant species among all the samples. Boxplot and violin plot will show every samples' relative
              abundance as a point.
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot.
            average_relative_abundance_threshold: (optional) Set a threshold, if you want to show all species with an
              equal or greater average relative abundance.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in the
              top. They will still be included in the relative abundances.
            prevalence_threshold: Prevalence threshold for a taxa to be kept in analysis.
            ascending: If set to True, from top to bottom, from least abundant taxa of the top to most abundant taxa.
        """
        plotting_options = kwargs.pop("plotting_options", {})
        mode = self._valid_mode_param(mode)
        if mode == "bargraph":
            fig = self._plot_most_abundant_taxa_bargraph(
                taxa_level=taxa_level,
                taxa_number=taxa_number,
                average_relative_abundance_threshold=average_relative_abundance_threshold,
                higher_classification=higher_classification,
                prevalence_threshold=prevalence_threshold,
                ascending=ascending,
                plotting_options=plotting_options,
                **kwargs,
            )
        else:
            fig = self._plot_most_what_taxa_boxplot_or_violin(
                "abundant",
                mode,
                taxa_level,
                taxa_number,
                determining_threshold=average_relative_abundance_threshold,
                higher_classification=higher_classification,
                threshold_on_other_variable=prevalence_threshold,
                ascending=ascending,
                plotting_options=plotting_options,
                **kwargs,
            )

        return fig

    def plot_sample_composition_most_abundant_taxa(
        self,
        taxa_level: str = "species",
        taxa_number: int = 20,
        average_relative_abundance_threshold: float = None,
        higher_classification: bool = True,
        prevalence_threshold: float = None,
        cluster_samples: bool = True,
        samples_order: List[str] = None,
        color_df: pd.DataFrame = None,
        sep_series: pd.Series = None,
        sep_how: str = None,
        **kwargs,
    ) -> go.Figure:
        """
        Plot taxa composition of samples for most abundant taxa.

        Args:
            taxa_level: Taxonomy level.
            taxa_number: Number of taxa to plot (skipped by average_relative_abundance_threshold).
            average_relative_abundance_threshold: (optional) Set a threshold, if you want to show all species with an
              equal or greater average relative abundance.
            higher_classification: Set to False, if you do not want OTU only defined at a higher level to appear in
              the top. They will still appear in "Others".
            prevalence_threshold: Prevalence threshold for a taxa to be kept in analysis.
            cluster_samples: Use clustering (skipped by samples_order).
            samples_order: List of samples to force ordering for visualization.
            color_df: Metadata to put as legend on the bottom of the graph.
            sep_series: Metadata used to order samples into subgroups (skipped by samples_order).
            sep_how: { None (default), 'color', 'labels' } Graphical way of showing the separation of the different
              subgroups (skipped if sep_series is empty/None).
        """
        data_df, taxa_number = self._compute_relative_abundances_taxa_dataframe(
            taxa_level=taxa_level,
            taxa_number=taxa_number,
            average_relative_abundance_threshold=average_relative_abundance_threshold,
            higher_classification=higher_classification,
            prevalence_threshold=prevalence_threshold,
            ascending=False,
        )[:2]
        data_df.loc["Others"] = 100 - data_df.sum()

        if data_df.shape[1] <= 1:  # only 1 sample, no need for ordering
            sep_series = None
        elif samples_order is not None:
            data_df = data_df.loc[:, samples_order]
        elif sep_series is not None:  # organize samples inside subgroups and concatenate subgroups one after another
            data_df, x_coor, subgps = self._divide_samples_into_subgroups_and_reorder(
                data_df, sep_series, cluster_samples=cluster_samples
            )
            if sep_how == "color":
                if color_df is None:
                    color_df = pd.DataFrame(sep_series)
                elif sep_series.name not in color_df.columns:
                    color_df = color_df.merge(pd.DataFrame(sep_series), right_index=True, left_index=True)
        elif cluster_samples:
            data_df = self._cluster_samples(data_df)

        # Make graph
        graph = MatrixBarGraph(data_df)
        # Plotting options
        title = f"{taxa_level.capitalize()} composition for the top {taxa_number} most abundant species across samples"
        if prevalence_threshold is not None:
            title += f" (present in at least {prevalence_threshold}% of samples)"

        default_plotting_options = {
            "layout": {
                "title": title,
                "xaxis_title": "Samples",
                "yaxis_title": "Relative abundance",
                "legend": {"traceorder": "normal"},
                "legend_title_text": "species",
            }
        }

        plotting_options = merge_dict(kwargs.pop("plotting_options", {}), default_plotting_options)

        if sep_series is not None and sep_how == "labels":
            show = kwargs.pop("show", True)
            output_file = kwargs.pop("output_file", False)
            if color_df is None:
                fig = graph.plot_one_graph(plotting_options=plotting_options, **kwargs, show=False)
            else:
                fig = graph.plot_complex_graph(color_df, plotting_options=plotting_options, **kwargs, show=False)
            fig = add_groups_annotations(fig, x_coor, subgps)
            graph._handle_output_plotly(fig, show, output_file)
        else:
            if color_df is None:
                fig = graph.plot_one_graph(plotting_options=plotting_options, **kwargs)
            else:
                fig = graph.plot_complex_graph(color_df, plotting_options=plotting_options, **kwargs)

        return fig
