import logging
from typing import List, Union

import pandas as pd

from moonstone.filtering.base import (
    BothAxisFiltering, CountsFiltering
)

logger = logging.getLogger(__name__)


class NoCountsFiltering(BothAxisFiltering, CountsFiltering):
    """
    Remove rows (default) or columns with no counts.
    """

    def filter(self) -> pd.DataFrame:
        indices = self.df.sum(axis=self.axis) != 0.0
        if self.axis == 1:
            return self.df.loc[indices]
        return self.df.loc[:, indices]


class NamesFiltering(BothAxisFiltering):
    """
    Filtering based on row (default) or column names.
    """

    def __init__(self, dataframe: pd.DataFrame, names: List[str], axis: int = 0, keep: bool = True):
        """
        :param names: list of row or column names
        :param axis: axis to apply filtering (index (0) or columns(1))
        :param keep: keep column (discard them if set to False)
        """
        self.names = list(names)
        self.keep = keep
        self._log_action()
        super().__init__(dataframe, axis=axis)

    def _log_action(self):
        if self.keep:
            logger.info(f"Selecting {self.names} from dataframe...")
        else:
            logger.info(f"Removing {self.names} from dataframe...")

    def _validate_parameters(self):
        if isinstance(self.df.index, pd.MultiIndex):
            error_message = f"{self.__class__.__name__} does not support filtering on MultiIndex dataframes." + \
                " You might want to use moonstone.filtering.TaxonomyNamesFiltering instead."
            raise TypeError(error_message)

    def _select_names(self):
        if self.axis == 0:
            return self.df.loc[self.names, :]
        return self.df.loc[:, self.names]

    def _exclude_names(self):
        return self.df.drop(self.names, axis=self.axis)

    def filter(self) -> pd.DataFrame:
        if self.keep:
            return self._select_names()
        else:
            return self._exclude_names()


class NaNPercentageFiltering(BothAxisFiltering):

    def __init__(self, dataframe: pd.DataFrame, percentage_nan_allowed: Union[int, float] = 80, axis: int = 0):
        """
        :param percentage_nan_allowed: maximum percentage of NaN values allowed (between 0 and 100)
        :param axis: axis to apply filtering (index (0) or columns(1))
        """
        self.percentage_nan_allowed = percentage_nan_allowed / 100
        super().__init__(dataframe, axis=axis)

    def filter(self) -> pd.DataFrame:
        if self.axis == 0:
            return self.df[self.df.isnull().mean(axis=1) <= self.percentage_nan_allowed]
        return self.df.loc[:, self.df.isnull().mean() <= self.percentage_nan_allowed]


class NumberOfDifferentValuesFiltering(BothAxisFiltering):
    def __init__(self, dataframe: pd.DataFrame, min_number_values: int = 2,
                 na: bool = False, axis: int = 0):
        """
        :param min_number_values: minimum number of different values accepted
        :param na: NaN values counted as a different value or not
        :param axis: axis to apply filtering (index (0) or columns(1))
        """
        self.min_number_values = min_number_values
        self.na = na
        super().__init__(dataframe, axis=axis)

    def filter(self) -> pd.DataFrame:
        new_df = self.df
        if self.axis == 0:
            for row in self.df.index:
                x = pd.Series(new_df.loc[row].unique())
                if not self.na:
                    x = x.dropna()
                if len(x) < self.min_number_values:
                    new_df.drop(row, inplace=True, axis=0)
            return new_df
        for col in self.df.columns:
            x = pd.Series(new_df[col].unique())
            if not self.na:
                x = x.dropna()
            if len(x) < self.min_number_values:
                new_df.drop(col, inplace=True, axis=1)
        return new_df
