
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

    def mean_gen_values(self):
        if getattr(self, "_mean_gen_values", None) is None:
            Log_Mean = self._log_df.mean(axis=1)
            setattr(self, "_mean_gen_values", Log_Mean)
        return self._mean_gen_values

    # Substracting the Log_mean to log original values
    def subs_log_values(self):
        if getattr(self, "_subs_log_values", None) is None:
            Substracted_df = self._log_df.sub(self._mean_gen_values, axis='rows')
            setattr(self, "_subs_log_values", Substracted_df)
        return self._subs_log_values

    # Calculating the median per sample
    def sample_log_median(self):
        if getattr(self, "_sample_log_median", None) is None:
            Median_df = self._subs_log_values.median()
            setattr(self, "_sample_log_median", Median_df)
        return self._sample_log_median

    def scaling_factor(self):
        if getattr(self, "_scaling_factor", None) is None:
            Scaling_factors = np.exp(self._sample_log_median)
            setattr(self, "_scaling_factor", Scaling_factors)
        return self._scaling_factor

    def scaled_df(self):
        if getattr(self, "_scaled_df", None) is None:
            final_df = self.df.div(self._scaling_factor)
            setattr(self, "_scaled_df", final_df)
        return self._scaled_df
