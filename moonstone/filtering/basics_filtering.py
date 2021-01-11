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
            logger.debug(f"Selecting {self.names} from dataframe...")
        else:
            logger.debug(f"Removing {self.names} from dataframe...")

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
        # intersection of names to remove or keep and index's/columns' names
        old_names = self.names
        if self.axis == 0:
            self.names = list(self.df.index.intersection(self.names))
        else:
            self.names = list(self.df.columns.intersection(self.names))
        if len(old_names) - len(self.names) > 0:
            list_names_not_found = list(set(old_names).difference(set(self.names)))
            list_names_not_found.sort()
            logger.info(f"{list_names_not_found}: {len(old_names) - len(self.names)} \
name(s) not found in the dataframe.")

        if self.keep:
            return self._select_names()
        else:
            return self._exclude_names()


class NaNPercentageFiltering(BothAxisFiltering):
    """
    Remove rows (default) or columns with a percentage of NaN values above a given percentage.
    """

    def __init__(self, dataframe: pd.DataFrame, percentage: Union[int, float] = 80, axis: int = 0):
        """
        :param percentage: maximum percentage of NaN values allowed (between 0 and 100)
        :param axis: axis to apply filtering (index (0) or columns(1))
        """
        self.percentage_of_nan_allowed = percentage
        super().__init__(dataframe, axis=axis)

    def filter(self) -> pd.DataFrame:
        thresh = self.df.shape[1-self.axis] - self.df.shape[1-self.axis] * (self.percentage_of_nan_allowed/100)
        return self.df.dropna(axis=self.axis, thresh=thresh)


class NumberOfDifferentValuesFiltering(BothAxisFiltering):
    """
    Filtering of rows (default) or columns based on the number of different (unique) values they hold.
    """

    def __init__(self, dataframe: pd.DataFrame,
                 min: int = None, max: int = None,
                 na: bool = False, axis: int = 0):
        """
        :param min: minimum number of different values accepted
        :param max: maximum number of different values accepted
        :param na: NaN values counted as a different value or not
        :param axis: axis to apply filtering (index (0) or columns(1))
        """
        if min is None and max is None:
            logger.warning("No min or max specified.")

        if min is None:
            self.min_number_values = 0
        else:
            self.min_number_values = min

        if max is None:
            self.max_number_values = float('inf')
        else:
            self.max_number_values = max

        self.na = na
        super().__init__(dataframe, axis=axis)

    def filter(self) -> pd.DataFrame:
        new_df = self.df
        if self.axis == 0:
            for row in self.df.index:
                x = pd.Series(new_df.loc[row].unique())
                if not self.na:
                    x = x.dropna()
                if len(x) < self.min_number_values or len(x) > self.max_number_values:
                    new_df.drop(row, inplace=True, axis=0)
            logger.info("%s/%s rows dropped", new_df.shape[0], self.df.shape[0])
            return new_df
        for col in self.df.columns:
            x = pd.Series(new_df[col].unique())
            if not self.na:
                x = x.dropna()
            if len(x) < self.min_number_values or len(x) > self.max_number_values:
                new_df.drop(col, inplace=True, axis=1)
        logger.info("%s/%s columns dropped", new_df.shape[1], self.df.shape[1])
        return new_df
