import logging
from typing import List

import pandas as pd

from moonstone.filtering.base import BothAxisFiltering

logger = logging.getLogger(__name__)


class NoCountsFiltering(BothAxisFiltering):
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
        if self.keep:
            return self._select_names()
        else:
            return self._exclude_names()
