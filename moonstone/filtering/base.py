from abc import ABC, abstractmethod

import pandas as pd

from moonstone.core.module_base import BaseModule, BaseDF


class BaseFiltering(BaseModule, BaseDF, ABC):

    def __init__(self, dataframe: pd.DataFrame):
        """
        :param dataframe: pandas dataframe, output of parsers' step
        """
        super().__init__(dataframe)
        self.raw_items_number = self.df.shape[0]
        self.raw_reads_number = self.df.sum().sum()

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
