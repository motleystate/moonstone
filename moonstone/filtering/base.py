import logging
from abc import ABC, abstractmethod

import pandas as pd

from moonstone.core.module_base import BaseModule, BaseDF

logger = logging.getLogger(__name__)


class BaseFiltering(BaseModule, BaseDF, ABC):

    @property
    def filtered_df(self) -> pd.DataFrame:
        """
        retrieves the filtered pandas dataframe
        """
        if getattr(self, '_filtered_df', None) is None:
            self._filtered_df = self.filter()
        return self._filtered_df

    @abstractmethod
    def filter(self) -> pd.DataFrame:
        """
        method that filters the items on your pandas dataframe.
        """
        pass


class CountsFiltering(BaseFiltering):

    def __init__(self, dataframe: pd.DataFrame):
        """
        :param dataframe: pandas dataframe, output of parsers' step
        """
        super().__init__(dataframe)
        self._check_counts()
        self.raw_items_number = self.df.shape[0]
        self.raw_reads_number = self.df.sum().sum()

    def _check_counts(self):
        """Stop on counts that are not numbers; warn about negative and blank counts. Nothing is changed."""
        non_numeric = [str(col) for col in self.df.columns if not pd.api.types.is_numeric_dtype(self.df[col])]
        if non_numeric:
            raise ValueError(f"Counts must be numbers, but these columns are not: {', '.join(non_numeric[:5])}"
                             f"{' ...' if len(non_numeric) > 5 else ''}. Check for text such as '12,5' or 'n/a'.")
        n_negative = int((self.df < 0).sum().sum())
        if n_negative:
            logger.warning("%i negative count(s) found: read counts are expected to be 0 or more.", n_negative)
        n_blank = int(self.df.isna().sum().sum())
        if n_blank:
            logger.warning("%i blank count cell(s) found. They are left blank; means use the other values only.",
                           n_blank)


class BothAxisFiltering(BaseFiltering, ABC):
    """
    Base for Filtering classes on both
    """

    def __init__(self, dataframe: pd.DataFrame, axis: int = 0):
        """
        :param axis: axis to apply filtering (index (0) or columns(1))
        """
        super().__init__(dataframe)
        self.axis = axis
        self._validate_parameters()

    def _validate_parameters(self):
        if self.axis not in [0, 1]:
            raise ValueError("axis needs to be 0 (index) or 1 (columns)")
