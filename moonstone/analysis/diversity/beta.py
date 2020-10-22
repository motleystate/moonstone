import logging
import re
from abc import ABC, abstractmethod

import pandas as pd
import skbio.diversity
from skbio.stats.distance import DistanceMatrix

from moonstone.analysis.diversity.base import DiversityBase

logger = logging.getLogger(__name__)


class BetaDiversity(DiversityBase, ABC):
    DIVERSITY_INDEXES_NAME = "beta_index"
    DEF_TITLE = "(beta diversity) distribution across the samples"

    def __init__(self, dataframe: pd.DataFrame):
        super().__init__(dataframe)
        self.index_name = " ".join(re.findall('[A-Z][^A-Z]*', self.__class__.__name__)).capitalize()

    @abstractmethod
    def compute_beta_diversity(self, df) -> DistanceMatrix:
        """
        method that compute the beta diversity
        """
        pass

    def compute_diversity(self) -> pd.Series:
        series = self.beta_diversity.to_series()
        series.name = self.DIVERSITY_INDEXES_NAME
        return series

    @property
    def beta_diversity(self):
        """
        DistanceMatrix from skbio.
        """
        if getattr(self, '_beta_diversity', None) is None:
            self._beta_diversity = self.compute_beta_diversity(self.df)
        return self._beta_diversity

    @property
    def beta_diversity_series(self):
        return self.diversity_indexes

    @property
    def beta_diversity_df(self):
        return self.beta_diversity.to_data_frame()

    def _get_grouped_df(self, metadata_series):
        df_list = []
        for group in metadata_series.unique():
            group_df = self.df.loc[:, metadata_series[metadata_series == group].index]
            beta_div_indexes = self.compute_beta_diversity(group_df).to_series().reset_index(drop=True).to_frame()
            beta_div_indexes.columns = [self.DIVERSITY_INDEXES_NAME]
            beta_div_indexes[metadata_series.name] = group
            df_list.append(beta_div_indexes)
        return pd.concat(df_list).dropna()


class BrayCurtis(BetaDiversity):
    """
    Perform calculation of the Bray Curtis for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df):    # compute_shannon_diversity
        """
        :param base: logarithm base chosen (NB : for ln, base=math.exp(1))
        """
        # steps to compute the index
        return skbio.diversity.beta_diversity("braycurtis", df.transpose(), df.columns)
