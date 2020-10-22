import logging
import re
from abc import ABC, abstractmethod
from typing import Union

import pandas as pd
import skbio

from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.plot.graphs.box import GroupBoxGraph, BoxGraph
from moonstone.plot.graphs.histogram import Histogram
from moonstone.plot.graphs.violin import GroupViolinGraph, ViolinGraph
from moonstone.utils.plot import (
    add_default_titles_to_plotting_options,
    add_x_to_plotting_options
)

logger = logging.getLogger(__name__)


class AlphaDiversity(BaseModule, BaseDF, ABC):
    ALPHA_DIVERSITY_INDEXES_NAME = "alpha_index"

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame]):
        """
        :param dataframe: taxonomy count table normalized to have each samples with the same number of reads
        """
        super().__init__(dataframe)
        self.index_name = " ".join(re.findall('[A-Z][^A-Z]*', self.__class__.__name__)).lower()

    @abstractmethod
    def compute_alpha_diversity(self) -> pd.Series:
        """
        method that compute the alpha diversity
        """
        pass

    @property
    def alpha_diversity_indexes(self):
        # call compute_alpha_diversity and store into self._alpha_indexes
        if getattr(self, '_alpha_diversity_indexes', None) is None:
            self._alpha_diversity_indexes = self.compute_alpha_diversity()
            self._alpha_diversity_indexes.name = self.ALPHA_DIVERSITY_INDEXES_NAME
        return self._alpha_diversity_indexes

    def _handle_plotting_options(self, plotting_options: dict, title: str, xlabel: str, ylabel: str, log_scale: bool):
        if plotting_options is None:
            plotting_options = {}
        if log_scale:
            plotting_options = add_x_to_plotting_options(plotting_options, 'yaxes', 'type', "log")
            ylabel = f"{ylabel} (log)"
        return add_default_titles_to_plotting_options(
            plotting_options, title, xlabel, ylabel
        )

    def _visualize_histogram(self, bins_size, plotting_options, show, output_file, log_scale: bool):
        title = self.index_name + " (alpha diversity) distribution across the samples"
        xlabel = self.index_name
        ylabel = "number of samples"
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        hist_fig = Histogram(self.alpha_diversity_indexes)
        hist_fig.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
        )

    def _visualize_violin(self, plotting_options: dict, show: bool, output_file: str, log_scale: bool):
        title = self.index_name + " (alpha diversity) distribution across the samples"
        xlabel = "Group"
        ylabel = self.index_name
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        violing_fig = ViolinGraph(self.alpha_diversity_indexes)
        violing_fig.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
        )

    def _visualize_boxplot(self, plotting_options: dict, show: bool, output_file: str, log_scale: bool):
        title = self.index_name + " (alpha diversity) distribution across the samples"
        xlabel = "Group"
        ylabel = self.index_name
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        violing_fig = BoxGraph(self.alpha_diversity_indexes)
        violing_fig.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
        )

    def visualize(
        self, mode: str = 'histogram', bins_size: Union[int, float] = 0.1, log_scale: bool = False,
        show: bool = True, output_file: str = False, plotting_options: dict = None,
    ):
        """
        :param mode: how to display (histogram, boxplot, or violin)
        :param bins_size: [mode histo only] size of the histo bins
        :param show: display your graph
        :param output_file: file path to output your html graph
        :param plotting_options: plotly plotting_options
        """
        if mode not in ['histogram', 'violin', 'boxplot']:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = "histogram"

        if mode == "histogram":
            self._visualize_histogram(bins_size, plotting_options, show, output_file, log_scale)
        elif mode == "violin":
            self._visualize_violin(plotting_options, show, output_file, log_scale)
        elif mode == "boxplot":
            self._visualize_boxplot(plotting_options, show, output_file, log_scale)

    def visualize_groups(
        self, metadata_df: pd.DataFrame, group_col: str,  mode: str = 'boxplot',
        log_scale: bool = False, colors: dict = None, groups: list = None,
        show: bool = True, output_file: str = False,
        plotting_options: dict = None,
    ):
        """
        :param metadata_df: dataframe containing metadata and information to group the data
        :param group_col: column from metadata_df used to group the data
        :param mode: how to display (boxplot, or violin)
        :param colors: overides color for groups. format {group_id: color}
        :param groups: specifically select groups to display among group_col
        :param show: display your graph
        :param output_file: file path to output your html graph
        :param plotting_options: plotly plotting_options
        """
        if mode not in ['violin', 'boxplot']:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = "boxplot"

        title = f"Distribution of <b>{self.index_name.capitalize()}</b> among samples<br><i>grouped by {group_col}"
        xlabel = f"{group_col}"
        ylabel = f"{self.index_name.capitalize()}"
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )

        df = pd.concat([metadata_df[group_col], self.alpha_diversity_indexes], axis=1).dropna()
        if mode == "violin":
            fig = GroupViolinGraph(df)
        elif mode == "boxplot":
            fig = GroupBoxGraph(df)
        fig.plot_one_graph(
            self.ALPHA_DIVERSITY_INDEXES_NAME, group_col,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            colors=colors,
            groups=groups,
        )


class ShannonIndex(AlphaDiversity):
    """
    Perform calculation of the shannon index for each samples of the dataframe
    """
    def compute_alpha_diversity(self, base: Union[int, float] = 2) -> pd.Series:    # compute_shannon_diversity
        """
        :param base: logarithm base chosen (NB : for ln, base=math.exp(1))
        """
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.shannon(self.df[i], base)
        return pd.Series(Seriesdic)


class SimpsonInverseIndex(AlphaDiversity):
    """
    Perform calculation of the simpson inverse index for each samples of the dataframe
    """
    def compute_alpha_diversity(self) -> pd.Series:
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.enspie(self.df[i])
        return pd.Series(Seriesdic)
