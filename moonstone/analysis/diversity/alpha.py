import logging
from abc import ABC
from typing import Union, Optional

import pandas as pd
import skbio

from moonstone.analysis.diversity.base import DiversityBase

from moonstone.analysis.statistical_test import (
    statistical_test_groups_comparison
)

logger = logging.getLogger(__name__)


class AlphaDiversity(DiversityBase, ABC):
    DIVERSITY_INDEXES_NAME = "alpha_index"
    DEF_TITLE = "(alpha diversity) distribution across the samples"

    @property
    def alpha_diversity_indexes(self):
        return self.diversity_indexes

    def compare_groups(
        self, metadata_df: pd.DataFrame, group_col: str, stat_test: str = 'mann_whitney_u',
        plotting_options: dict = None,
        show_visualization: Optional[bool] = False, output_visualization_file: Optional[str] = False
    ):
        self.stat_test_group_matrix = statistical_test_groups_comparison(
            self.alpha_diversity_indexes, metadata_df[group_col], stat_test=stat_test
            )
        if show_visualization or output_visualization_file:
            self.visualize_groups(metadata_df, group_col, plotting_options=plotting_options,
                                  show=show_visualization, output_file=output_visualization_file)
        return self.stat_test_group_matrix


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
