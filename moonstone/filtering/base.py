import pandas as pd

from moonstone.core.module_base import BaseModule


class BaseFiltering(BaseModule):

    def __init__(self, dataframe: pd.DataFrame):
        """
        :param dataframe: pandas dataframe, output of parsers' step
        """
        self.df = dataframe
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

    def filter(self) -> pd.DataFrame:
        """
        method that filter the items on your pandas dataframe
        """
        return self.df
