from abc import ABC

from typing import Union, Optional

import pandas as pd
import re
import skbio

from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.plot.graphs.histogram import Histogram
from moonstone.utils.plot import (
    add_default_titles_to_plotting_options
)


class BaseAnalysis(BaseModule, BaseDF, ABC):
    """
    ????
    """
    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame]):
        super().__init__(dataframe)


class AlphaDiversity(BaseModule, BaseDF, ABC):

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame]):
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
        return self._alpha_diversity_indexes

    def visualize(self, bins_size: Union[int, float] = 0.1, plotting_options: dict = None,
                  show: Optional[bool] = True, output_file: Optional[str] = False):

        title = self.index_name+" (alpha diversity) distribution across the samples"
        xlabel = self.index_name
        ylabel = "number of samples"

        if plotting_options is None:
            plotting_options = {'layout': {'title_text': title, 'title_x': 0.5},
                                'xaxes': {'title_text': xlabel},
                                'yaxes': {'title_text': ylabel}}
        else:
            plotting_options = add_default_titles_to_plotting_options(plotting_options,
                                                                      title,
                                                                      xlabel, ylabel)

        hist_fig = Histogram(self.alpha_diversity_indexes)
        hist_fig.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            )


class ShannonIndex(AlphaDiversity):

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

    def compute_alpha_diversity(self):
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.enspie(self.df[i])
        return pd.Series(Seriesdic)
