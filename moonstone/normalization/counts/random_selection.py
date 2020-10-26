import logging

import numpy as np
import pandas as pd

from moonstone.filtering.basics_filtering import NamesFiltering
from moonstone.normalization.base import BaseNormalization

logger = logging.getLogger(__name__)


class RandomSelection(BaseNormalization):
    """
    Randomly select a given number of counts (threshold) among all different items
    (genes, taxonomical annotation...) for each sample. Random selection takes into account the
    initial counts to influence the probability of picking one item or another.
    """

    def __init__(self, df: pd.DataFrame, threshold: int = None, random_seed: int = 2935):
        """
        :param threshold: total number of counts to pick by sample
        :param random_seed: random seed to use for random picking of counts
        """
        super().__init__(df)
        self.random_seed = random_seed
        if threshold is not None:
            self.threshold = threshold
        else:
            self.threshold = int(self.df.sum().min())
        # Filters out samples below the threshold
        self.samples_to_remove = self.raw_df.columns[self.raw_df.sum() < self.threshold]
        if not self.samples_to_remove.empty:
            self.df = NamesFiltering(self.raw_df, self.samples_to_remove, axis=1, keep=False).filtered_df

    def _randomly_select_counts(self, column_name: str):
        np.random.seed(self.random_seed)  # set the random seed
        counts = self.df[column_name]
        if counts.sum() <= self.threshold + 1:
            return counts
        probabilities = counts / counts.sum()
        new_counts = np.unique(
            np.random.choice(counts.index, self.threshold, p=probabilities), return_counts=True
        )
        return pd.Series(new_counts[1], index=new_counts[0])

    def normalize(self) -> pd.DataFrame:
        normalized_df = pd.DataFrame()
        cpt = 0
        total = len(self.df.columns)
        for sample in self.df.columns:
            normalized_df = pd.concat([normalized_df, self._randomly_select_counts(sample)], axis=1)
            cpt += 1
            if cpt % 10 == 0:
                logger.info(f"{cpt}/{total} done so far...")
        logger.info(f"[Done] {cpt}/{total}.")
        normalized_df.columns = self.df.columns
        return normalized_df.fillna(0).astype(float)


class TaxonomyRandomSelection(RandomSelection):
    """
    Allow random selection for taxonomy multi-indexed dataframes.
    """

    def __init__(self, df: pd.DataFrame, concat_char: str = ';', *args, **kwargs):
        self.concat_char = concat_char
        no_index_df = df.reset_index()
        self.index_names = df.index.names
        new_df = no_index_df.set_index(
            no_index_df[self.index_names].agg(self.concat_char.join, axis=1)
        ).drop(self.index_names, axis=1)
        super().__init__(new_df, **kwargs)

    def normalize(self) -> pd.DataFrame:
        single_index_norm_df = super().normalize()
        multi_index_norm_df = single_index_norm_df.reset_index(drop=True)
        multi_index_norm_df.index = single_index_norm_df.index.str.split((self.concat_char), expand=True)
        multi_index_norm_df.index.names = self.index_names
        return multi_index_norm_df
