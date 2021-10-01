import logging
from abc import ABC
from typing import Union

import pandas as pd
import skbio

from moonstone.analysis.diversity.base import (
    DiversityBase, PhylogeneticDiversityBase
)

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
    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame], base: Union[int, float] = 2):
        """
        :param base: logarithm base chosen (NB : for ln, base=math.exp(1))
        """
        super().__init__(dataframe)
        self.base = base

    def compute_diversity(self) -> pd.Series:    # compute_shannon_diversity
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.shannon(self.df[i], base=self.base)
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


class Chao1Index(AlphaDiversity):
    """
    Perform calculation of the Choa1 index for each samples of the dataframe
    """
    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame], bias_corrected: bool = True):
        """
        :param bias_corrected
        """
        super().__init__(dataframe)
        self.bias_corrected = bias_corrected

    def compute_diversity(self) -> pd.Series:
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.chao1(self.df[i], bias_corrected=self.bias_corrected)
        return pd.Series(Seriesdic)


class FaithsPhylogeneticDiversity(AlphaDiversity, PhylogeneticDiversityBase):
    """
    Perform calculation of the Faith's PD for each samples of the dataframe
    """
    def compute_diversity(self) -> pd.Series:
        # steps to compute the index
        seriesdic = {}
        otu_ids = self.df.index

        missing_ids = self._verification_otu_ids_in_tree(otu_ids)
        if len(missing_ids) > 0:
            if not self.force_computation:
                raise RuntimeError(f"INCOMPLETE TREE: missing {missing_ids}.")
            else:
                logger.warning(f"INCOMPLETE TREE: missing {missing_ids}.\n\
Computation of the Faith's diversity using only the OTU IDs present in the Tree.")
                otu_ids = list(set(otu_ids) - set(missing_ids))

        for i in self.df.columns:
            seriesdic[i] = skbio.diversity.alpha.faith_pd(
                self.df[i].loc[otu_ids], otu_ids, self.tree, validate=self.validate
                )
        return pd.Series(seriesdic)
