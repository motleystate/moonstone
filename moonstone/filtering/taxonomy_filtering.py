from typing import List

import pandas as pd

from moonstone.filtering.base import CountsFiltering


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
    taxonomy below this mean.
    """

    def __init__(
        self, dataframe: pd.DataFrame, mean_value: float,
        level: str = 'species'
    ):
        """
        :param mean_value: mean among all samples to be kept
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
        level_df = self.df.groupby(self.level).sum()
        selected_names = level_df[level_df.mean(axis=1) > self.mean_value].index
        corresponding_rows = self.df.index.get_level_values(self.level).isin(selected_names)
        return self.df.loc[corresponding_rows]
