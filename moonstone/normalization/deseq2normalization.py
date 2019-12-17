import numpy as np
import math
import logging
import pandas as pd


class DESeq2Normalization:

    def __init__(self, df, log_number=np.e, zero_threshold=80, normalization_level=0):
        self.df = df
        self.log = log_number
        self.zero_threshold = zero_threshold
        self.normalization_level = normalization_level
        if isinstance(df.index, pd.core.index.MultiIndex):
            self.grouped_df = df.groupby(level=self.normalization_level).sum()
        else:
            self.grouped_df = df

    def non_zero_df(self, df):
        """
        This function removes rows with 0 reads
        """
        threshold = math.ceil(df.shape[1] * self.zero_threshold/100)
        non_zero_dataf = df.replace(0, np.nan).dropna(thresh=threshold).astype('float')
        total_len = len(df)
        non_zero_df_len = len(non_zero_dataf)
        if non_zero_df_len / total_len * 100 <= 50:
            logging.warning("{} rows were dropped, which represents {} % of the sample".format(
                            total_len - non_zero_df_len, (total_len - non_zero_df_len) / total_len*100))
        self._removed_zero_df = df[~df.index.isin(non_zero_dataf.index)]
        return non_zero_dataf

    def log_df(self, df):
        return df.applymap(lambda x: math.log(x, self.log))

    @property
    def removed_zero_df(self):
        """
        gives the data frame with the rows that were removed for having too many zeros.
        this attribute is computed during the non_zero_df function
        """
        if getattr(self, "_removed_zero_df", None) is None:
            logging.warning("Computing the scaling factors beforehand is required")
            return None
        else:
            return self._removed_zero_df

    def remove_zero_and_log(self, df):
        return self.log_df(self.non_zero_df(df))

    def calculating_and_substracting_mean_row(self, df):
        """
        Substracting the mean row to original values
        """
        return df.sub(df.mean(axis=1), axis='rows')

    @property
    def scaling_factors(self):
        if getattr(self, "_scaling_factors", None) is None:
            non_zero_log_df = self.remove_zero_and_log(self.grouped_df)
            substracted_mean_df = self.calculating_and_substracting_mean_row(non_zero_log_df)
            scaling_factors = substracted_mean_df.rpow(self.log).replace(np.nan, 0).median()
            setattr(self, "_scaling_factors", scaling_factors)
        return self._scaling_factors

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            final_df = self.df.div(self.scaling_factors)
            setattr(self, "_normalized_df", final_df)
        return self._normalized_df
