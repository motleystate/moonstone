import numpy as np
import math


class DESeq2Normalization:

    def __init__(self, df, log_number=np.e):
        self.df = df
        self.log = log_number

    def non_zero_df(self, df):
        """
        This function removes rows with 0 reads
        """
        return df.replace(0, np.nan).dropna().astype(int)

    def log_df(self, df):
        return df.applymap(lambda x: math.log(x, self.log))

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
            non_zero_log_df = self.remove_zero_and_log(self.df)
            substracted_mean_df = self.calculating_and_substracting_mean_row(non_zero_log_df)
            scaling_factors = substracted_mean_df.applymap(lambda x: math.pow(self.log, x)).median()
            setattr(self, "_scaling_factors", scaling_factors)
        return self._scaling_factors

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            final_df = self.df.div(self.scaling_factors)
            setattr(self, "_normalized_df", final_df)
        return self._normalized_df
