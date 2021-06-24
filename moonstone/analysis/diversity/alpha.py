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


class Chao1Index(AlphaDiversity):
    """
    Perform calculation of the Choa1 index for each samples of the dataframe
    """
    def compute_diversity(self, bias_corrected: bool = True) -> pd.Series:
        # steps to compute the index
        Seriesdic = {}
        for i in self.df.columns:
            Seriesdic[i] = skbio.diversity.alpha.chao1(self.df[i], bias_corrected)
        return pd.Series(Seriesdic)


class FaithsPhylogeneticDiversity(AlphaDiversity):
    """
    Perform calculation of the Faith's PD for each samples of the dataframe
    """
    def __init__(
        self,
        taxonomy_dataframe: pd.DataFrame,
        taxonomy_tree: skbio.TreeNode
    ):
        super().__init__(taxonomy_dataframe)
        self.tree = taxonomy_tree

    def compute_diversity(self, validate: bool = True) -> pd.Series:
        """
        Args:
            validate: skbio argument. "If False, validation of the input won’t be performed.
            This step can be slow, so if validation is run elsewhere it can be disabled here.
            However, invalid input data can lead to invalid results or error messages that
            are hard to interpret, so this step should not be bypassed if you’re not certain
            that your input data are valid. See skbio.diversity for the description of what
            validation entails so you can determine if you can safely disable validation.
        """
        # steps to compute the index
        seriesdic = {}
        otu_ids = self.df.index

        missing_ids = []
        for otu_id in otu_ids:
            try:
                self.tree.find(otu_id)
            except:
                missing_ids += [otu_id]
        if len(missing_ids) > 0:
            raise RuntimeError(f"INCOMPLETE TREE: missing {missing_ids}.")

        for i in self.df.columns:
            seriesdic[i] = skbio.diversity.alpha.faith_pd(self.df[i], otu_ids, self.tree, validate=validate)
        return pd.Series(seriesdic)
