import logging

import numpy as np
import pandas as pd

from moonstone.normalization.base import BaseNormalization

logger = logging.getLogger(__name__)


class RandomSelection(BaseNormalization):
    """
    Randomly select counts for each sample based on a threshold.
    """

    def __init__(self, df: pd.DataFrame, threshold: int = None, random_seed: int = 2935):
        super().__init__(df)
        self.random_seed = random_seed
        if threshold is not None:
            self.threshold = threshold
        else:
            self.threshold = self.df.sum().min()
        # Filters out samples below this threshold?

    def _randomly_select_counts(self, column_name: str):
        np.random.seed(self.random_seed)  # set the random seed
        counts = self.raw_df[column_name]
        if counts.sum() <= self.threshold:
            return counts
        probabilities = counts / counts.sum()
        new_counts = np.unique(
            np.random.choice(counts.index, self.threshold, p=probabilities), return_counts=True
        )
        return pd.Series(new_counts[1], index=new_counts[0])

    def normalize(self) -> pd.DataFrame:
        normalized_df = pd.DataFrame()
        for sample in self.raw_df.columns:
            normalized_df = pd.concat([normalized_df, self._randomly_select_counts(sample)], axis=1)
        normalized_df.columns = self.raw_df.columns
        return normalized_df.fillna(0).astype(int)
