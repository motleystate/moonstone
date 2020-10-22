from abc import ABC, abstractmethod

from skbio.diversity import beta_diversity
from skbio.stats.distance import DistanceMatrix

from moonstone.core.module_base import BaseModule, BaseDF


class BetaDiversity(BaseModule, BaseDF, ABC):
    BETA_DIVERSITY_INDEXES_NAME = "beta_index"

    @abstractmethod
    def compute_beta_diversity(self) -> DistanceMatrix:
        """
        method that compute the beta diversity
        """
        pass

    @property
    def beta_diversity(self):
        """
        DistanceMatrix from skbio.
        """
        if getattr(self, '_beta_diversity', None) is None:
            self._beta_diversity = self.compute_beta_diversity(self.df)
            self._beta_diversity.name = self.BETA_DIVERSITY_INDEXES_NAME
        return self._beta_diversity

    @property
    def beta_diversity_series(self):
        return self.beta_diversity.to_series()

    @property
    def beta_diversity_df(self):
        return self.beta_diversity.to_data_frame()


class BrayCurtis(BetaDiversity):
    """
    Perform calculation of the bray curtis for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df):    # compute_shannon_diversity
        """
        :param base: logarithm base chosen (NB : for ln, base=math.exp(1))
        """
        # steps to compute the index
        return beta_diversity("braycurtis", df.transpose(), df.columns)
