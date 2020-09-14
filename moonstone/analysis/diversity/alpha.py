from abc import ABC

from typing import Union, Optional

import pandas as pd
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


class ShannonIndex(BaseModule, BaseDF, ABC):

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame]):
        super().__init__(dataframe)

    def compute_shannon_diversity(self, base: Union[int, float] = 2):    # compute_shannon_diversity
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.shannon(self.df[i], base)
        return pd.Series(Seriesdic)

    @property
    def shannon_indexes(self):
        # call compute_shannon_diversity and store into self._shannon_indexes
        if getattr(self, '_shannon_indexes', None) is None:
            self._shannon_indexes = self.compute_shannon_diversity()
        return self._shannon_indexes

    def visualize(self, bins_size: Union[int, float] = 0.1, plotting_options: dict = None,
                  show: Optional[bool] = True, output_file: Optional[str] = False):

        title = "shannon indexes (alpha diversity) distribution across the samples"
        xlabel = "shannon index"
        ylabel = "number of samples"

        if plotting_options is None:
            plotting_options = {'layout': {'title_text': title, 'title_x': 0.5},
                                'xaxes': {'title_text': xlabel},
                                'yaxes': {'title_text': ylabel}}
        else:
            plotting_options = add_default_titles_to_plotting_options(plotting_options,
                                                                      title,
                                                                      xlabel, ylabel)

        hist_fig = Histogram(self._shannon_indexes)
        hist_fig.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            )
