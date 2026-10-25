import logging
from typing import List

import pandas as pd

from moonstone.filtering.base import CountsFiltering

logger = logging.getLogger(__name__)


class TaxonomyNamesFiltering(CountsFiltering):
    """
    Filtering a Taxonomy multiindexed dataframe on index names at a chosen level.
    """

    def __init__(
        self, dataframe: pd.DataFrame, names: List[str],
        level: str = 'species', keep: bool = True
    ):
        """
        :param names: list of index names
        :param level: level of the MultiIndex to filter on
        :param keep: keep column (discard them if set to False)
        """
        super().__init__(dataframe)
        self.names = names
        self.keep = keep
        self.level = level
        self._validate_parameters()

    def _validate_parameters(self):
        if self.level not in self.df.index.names:
            error_message = f"{self.level} not a valid level. Must be among {self.df.index.names}"
            raise ValueError(error_message)

    def filter(self) -> pd.DataFrame:
        corresponding_rows = self.df.index.get_level_values(self.level).isin(self.names)
        if self.keep:
            return self.df.loc[corresponding_rows]
        else:
            return self.df.loc[~corresponding_rows]


class TaxonomyMeanFiltering(CountsFiltering):
    """
    Filtering a Taxonomy multiindexed dataframe on sample mean at a chosen level.

    This means you select a mean value for all samples and it will discard all selected
    taxonomy below this mean. Taxa whose mean equals the mean value are kept.
    Rows with no name at the chosen level are summed together as one unnamed group, not dropped.
    """

    def __init__(
        self, dataframe: pd.DataFrame, mean_value: float,
        level: str = 'species'
    ):
        """
        :param mean_value: minimum mean among all samples for a taxon to be kept
        :param level: level of the MultiIndex to filter on
        """
        super().__init__(dataframe)
        self.mean_value = mean_value
        self.level = level
        self._validate_parameters()

    def _validate_parameters(self):
        if self.level not in self.df.index.names:
            error_message = f"{self.level} not a valid level. Must be among {self.df.index.names}"
            raise ValueError(error_message)

    def filter(self) -> pd.DataFrame:
        # Rows with no name at this level are kept together as one blank-named group, not dropped.
        n_unnamed = int(self.df.index.get_level_values(self.level).isna().sum())
        if n_unnamed:
            logger.warning("%i row(s) have no %s name; they are filtered together as one unnamed group.",
                           n_unnamed, self.level)
        level_df = self.df.groupby(level=self.level, dropna=False).sum()
        selected_names = level_df[level_df.mean(axis=1) >= self.mean_value].index
        corresponding_rows = self.df.index.get_level_values(self.level).isin(selected_names)
        return self.df.loc[corresponding_rows]
