import logging
import re
from abc import ABC, abstractmethod
from typing import Union, Optional

import pandas as pd
import skbio

from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.plot.graphs.box import GroupBoxGraph, BoxGraph
from moonstone.plot.graphs.histogram import Histogram
from moonstone.plot.graphs.violin import GroupViolinGraph, ViolinGraph
from moonstone.utils.plot import (
    add_default_titles_to_plotting_options
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
    def compute_alpha_diversity(self):
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

    def _visualize_histogram(self, bins_size, plotting_options, show, output_file):
        title = self.index_name + " (alpha diversity) distribution across the samples"
        xlabel = self.index_name
        ylabel = "number of samples"
        plotting_options = add_default_titles_to_plotting_options(
            plotting_options, title, xlabel, ylabel
        )
        hist_fig = Histogram(self.alpha_diversity_indexes)
        hist_fig.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
        )

    def _visualize_violin(self, plotting_options, show, output_file, log_scale: bool = False):
        title = self.index_name + " (alpha diversity) distribution across the samples"
        xlabel = "Group"
        ylabel = self.index_name
        plotting_options = add_default_titles_to_plotting_options(
            plotting_options, title, xlabel, ylabel
        )
        violing_fig = ViolinGraph(self.alpha_diversity_indexes)
        violing_fig.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            log_scale=log_scale,
        )

    def _visualize_boxplot(self, plotting_options, show, output_file, log_scale: bool = False):
        title = self.index_name + " (alpha diversity) distribution across the samples"
        xlabel = "Group"
        ylabel = self.index_name
        plotting_options = add_default_titles_to_plotting_options(
            plotting_options, title, xlabel, ylabel
        )
        violing_fig = BoxGraph(self.alpha_diversity_indexes)
        violing_fig.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            log_scale=log_scale,
        )

    def visualize(self, mode: str = 'histogram', bins_size: Union[int, float] = 0.1, plotting_options: dict = None,
                  show: Optional[bool] = True, output_file: Optional[str] = False, log_scale: bool = False):
        """
        :param mode: how to display (histogram or violin)
        :param bins_size: [mode histo only] size of the histo bins
        """
        if mode not in ['histogram', 'violin', 'boxplot']:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = "histogram"

        if mode == "histogram":
            self._visualize_histogram(bins_size, plotting_options, show, output_file)
        elif mode == "violin":
            self._visualize_violin(plotting_options, show, output_file, log_scale=log_scale)
        elif mode == "boxplot":
            self._visualize_boxplot(plotting_options, show, output_file, log_scale=log_scale)

    def visualize_groups(
        self, metadata_df: pd.DataFrame, group_col: str, plotting_options: dict = None,
        show: Optional[bool] = True, output_file: Optional[str] = False,
        log_scale: bool = False, mode: str = 'boxplot', colors: dict = None,
    ):
        """
        :param metadata_df: dataframe containing metadata and information to group the data
        :param group_col: column from metadata_df used to group the data
        """
        if mode not in ['violin', 'boxplot']:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = "boxplot"

        title = f"Distribution of <b>{self.index_name.capitalize()}</b> among samples<br><i>grouped by {group_col}"
        xlabel = f"{group_col}"
        ylabel = f"{self.index_name.capitalize()}"
        plotting_options = add_default_titles_to_plotting_options(
            plotting_options, title, xlabel, ylabel
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
            log_scale=log_scale,
            colors=colors,
        )


class ShannonIndex(AlphaDiversity):
    """
    Perform calculation of the shannon index for each samples of the dataframe
    """
    def compute_alpha_diversity(self, base: Union[int, float] = 2):    # compute_shannon_diversity
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
    def compute_alpha_diversity(self):
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.enspie(self.df[i])
        return pd.Series(Seriesdic)
