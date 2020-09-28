from abc import ABC

from typing import Union, Optional

import pandas as pd
import re
from scipy import stats
import skbio

from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.plot.graphs.histogram import Histogram
from moonstone.utils.plot import (
    add_default_titles_to_plotting_options
)


class AlphaDiversity(BaseModule, BaseDF, ABC):

    bins_size = 0.1

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
        return self._alpha_diversity_indexes

    def stats(self):
        # copy-paste from Descriptive in moonstone.analysis.stats
        print(f"Descriptive statistics of {self.index_name}:")

        properties = ['Number of Variable', 'Min / Max', 'Mean', 'Variance', 'Skewness', 'Kurtosis / Fisher']
        for i in range(len(stats.describe(self.alpha_diversity_indexes))):
            print(f"\t{properties[i]} = {stats.describe(self.alpha_diversity_indexes)[i]}")

    def visualize(self, bins_size: Union[int, float] = None, plotting_options: dict = None,
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

        if bins_size is None:
            bins_size = self.bins_size

        hist_fig = Histogram(self.alpha_diversity_indexes)
        hist_fig.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            )


class ShannonIndex(AlphaDiversity):
    """
    Perform calculation of the shannon index for each samples of the dataframe
    """
    bins_size = 0.1

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
    bins_size = 100

    def compute_alpha_diversity(self):
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.enspie(self.df[i])
        return pd.Series(Seriesdic)
