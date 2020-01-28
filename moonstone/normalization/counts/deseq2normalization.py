import numpy as np
import math
import logging
import pandas as pd

from moonstone.normalization.base import BaseNormalization


class DESeq2Normalization(BaseNormalization):
    """
    Brief explanation of how this normalization module works:
    "gen" "spec" 1  2  3                  "gen"  1    2     3                 "gen"    "spec"   1       2        3
    gen_1  spec_1 a  b  c  groupby('genus') gen_1 a+d  b+f  c+g   normalization gen_1  spec_1 norm(a)  norm(b)  norm(c)
    gen_1  spec_2 d  f  g    .sum()                                 --------->  gen_1  spec_2 norm(d)  norm(f)  norm(g)
    gen_2  spec_3 j  m  o  ---------------> gen_2 j+k  m+l  o+p                 gen_2  spec_3 norm(j)  norm(m)  norm(o)
    gen_2  spec_4 k  l  p                                                       gen_2  spec_4 norm(k)  norm(l)  norm(p)
                                            scaling factors based
                                            on this table, but applied
                                            to the original df
"""

    def __init__(self, df, log_number=np.e, zero_threshold=80, normalization_level=0):
        """
        :param normalization_level: At which level of a multi-index you want the normalization to be perfomed
        """
        super().__init__(df)
        self.log = log_number
        self.zero_threshold = zero_threshold
        self.normalization_level = normalization_level
        if isinstance(df.index, pd.core.index.MultiIndex):
            self.grouped_df = df.groupby(level=self.normalization_level).sum()
        else:
            self.grouped_df = df

    def non_zero_df(self, df):
        """
        This method removes rows with 0 reads
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
        gives the dataframe with the rows that were removed for having too many zeros.
        this attribute is computed during the non_zero_df function
        """
        if getattr(self, "_removed_zero_df", None) is None:
            logging.warning("Computing the scaling factors beforehand is required to access this dataframe")
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

    def normalize(self):
        return self.df.div(self.scaling_factors)
