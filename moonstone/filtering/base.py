from moonstone.core.module_base import BaseModule


class BaseFiltering(BaseModule):

    def __init__(self, dataframe):
        self.df = dataframe
        self.steps = []
        self.raw_items_number = self.df.shape[0]
        self.raw_reads_number = self.df.sum().sum()

    @property
    def filtered_df(self):
        if getattr(self, '_filtered_df', None) is None:
            self._filtered_df = self.filter()
        return self._filtered_df

    def filter(self):
        return self.df
