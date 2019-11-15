import numpy as np


class DESeq2_normalization:

    def __init__(self, df):
        self.df = df

    def non_zero_df(self, df):
        """
        This function removes rows with 0 reads
        """
        return df.replace(0, np.nan).dropna().astype(int)

    def log_df(self, df, log_nb=None):
        if log_nb is None:
            return np.log(df)
        else:
            return np.log(df) / np.log(log_nb)

    def remove_zero_and_log(self, df):
        return self.log_df(self.non_zero_df(df))

    def calculating_and_substracting_mean_row(self, df):
        """
        Substracting the mean row to original values
        """
        return df.sub(df.mean(axis=1), axis='rows')

    # Only missing how to change the base_exponent
    @property
    def compute_scaling_factor(self):
        if getattr(self, "_compute_scaling_factor", None) is None:
            non_zero_log_df = self.remove_zero_and_log(self.df)
            substracted_mean_df = self.calculating_and_substracting_mean_row(non_zero_log_df)
            Scaling_factors = np.exp(substracted_mean_df.median())
            setattr(self, "_compute_scaling_factor", Scaling_factors)
        return self._compute_scaling_factor

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            final_df = self.df.div(self.compute_scaling_factor)
            setattr(self, "_normalized_df", final_df)
        return self._normalized_df
