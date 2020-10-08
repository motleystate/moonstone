import logging
import re
from abc import ABC
from typing import Union, Optional

import pandas as pd
import skbio
from plotly import graph_objects as go

from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.plot.graphs.histogram import Histogram
from moonstone.plot.graphs.violin import ViolinGraph
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

    def compute_alpha_diversity(self):
        """
        method that compute the alpha diversity
        """
        return "choose index computation method"

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

    def _visualize_violin(self, plotting_options, show, output_file):
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
        )

    def visualize(self, mode: str = 'histogram', bins_size: Union[int, float] = 0.1, plotting_options: dict = None,
                  show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        :param mode: how to display (histogram or violin)
        :param bins_size: [mode histo only] size of the histo bins
        """
        if mode not in ['histogram', 'violin']:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = "histogram"

        if mode == "histogram":
            self._visualize_histogram(bins_size, plotting_options, show, output_file)
        elif mode == "violin":
            self._visualize_violin(plotting_options, show, output_file)

    def visualize_groups(self, metadata_df: pd.DataFrame, group_col: str, plotting_options: dict = None,
                         show: Optional[bool] = True, output_file: Optional[str] = False):
        groups = list(metadata_df[group_col].unique())
        groups.sort()

        df = pd.concat([metadata_df['MOTIF_INC'], self.alpha_diversity_indexes], axis=1).dropna()

        fig = go.Figure()
        for group in groups:
            fig.add_trace(go.Violin(x=df['MOTIF_INC'][df['MOTIF_INC'] == group],
                                    y=df['alpha_index'][df['MOTIF_INC'] == group],
                                    name=str(group),
                                    box_visible=True,
                                    points='all',
                                    meanline_visible=True,
                                    text=df.index))
        fig.update_layout(
            title_text=f"<b>{self.index_name.capitalize()}-index</b><br><i>grouped by {group_col}",
            xaxis_title=f"{group_col}",
            yaxis_title="Alpha indexes"
        )
        fig.show()



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
