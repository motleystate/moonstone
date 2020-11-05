import logging
from abc import ABC
from typing import Union

import pandas as pd
import skbio

from moonstone.analysis.diversity.base import DiversityBase

logger = logging.getLogger(__name__)


class AlphaDiversity(DiversityBase, ABC):
    DIVERSITY_INDEXES_NAME = "alpha_index"
    DEF_TITLE = "(alpha diversity) distribution across the samples"

    @property
    def alpha_diversity_indexes(self):
        return self.diversity_indexes


class ShannonIndex(AlphaDiversity):
    """
    Perform calculation of the shannon index for each samples of the dataframe
    """
    def compute_diversity(self, base: Union[int, float] = 2) -> pd.Series:    # compute_shannon_diversity
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
    def compute_diversity(self) -> pd.Series:
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.enspie(self.df[i])
        return pd.Series(Seriesdic)
