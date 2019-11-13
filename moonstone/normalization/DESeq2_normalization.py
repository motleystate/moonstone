import numpy as np


class DESeq2_normalization:

    def __init__(self, df):
        self.df = df

    def non_zero_df(self):
        """
        This function removes lines with 0 reads
        """
        return self.df.copy().replace(0, np.nan).dropna().astype(int)

    def log_df(self):
        return np.log(self.non_zero_df())

    def mean_gen_values(self):
        return self.log_df().mean(axis=1)

    def subs_log_values(self):
        """
        Substracting the Log_mean to log original values
        """
        return self.log_df().sub(self.mean_gen_values(), axis='rows')

    def sample_log_median(self):
        """
        Calculating the median per sample
        """
        return self.subs_log_values().median()

    @property
    def scaling_factor(self):
        if getattr(self, "_scaling_factor", None) is None:
            Scaling_factors = np.exp(self.sample_log_median())
            setattr(self, "_scaling_factor", Scaling_factors)
        return self._scaling_factor

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            final_df = self.df.div(self.scaling_factor)
            setattr(self, "_normalized_df", final_df)
        return self._normalized_df
