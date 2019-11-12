import pandas as pd
import numpy as np


class DESeq2_normalization:

    def __init__(self, df):
        self.df = df

    def non_zero_df(self):
        if getattr(self, "_non_zero_df", None) is None:
            df_without_zero = self.df.copy().replace(0, np.nan).dropna()
            setattr(self, "_non_zero_df", df_without_zero)
        return self._non_zero_df

    def log_df(self):
        if getattr(self, "_log_df", None) is None:
            log_non_zero_df = np.log(self._non_zero_df)
            setattr(self, "_log_df", log_non_zero_df)
        return self._log_df
